import argparse
import math
import os
import sys
import pandas as pd

# trunk-ignore(flake8/E402)
import tensorflow.compat.v1 as tf

tf.disable_v2_behavior()


class ImageReader:
    """Helper class that provides TensorFlow image coding utilities."""

    def __init__(self):
        # Initializes function that decodes RGB JPEG data.
        self._decode_jpeg_data = tf.placeholder(dtype=tf.string)
        self._decode_jpeg = tf.image.decode_jpeg(self._decode_jpeg_data, channels=3)

    def read_image_dims(self, sess, image_data):
        image = self.decode_jpeg(sess, image_data)
        return image.shape[0], image.shape[1]

    def decode_jpeg(self, sess, image_data):
        image = sess.run(
            self._decode_jpeg, feed_dict={self._decode_jpeg_data: image_data}
        )
        assert len(image.shape) == 3
        assert image.shape[2] == 3
        return image


def _get_dataset_filename(output_dir, split_name, shard_id, NUM_SHARDS):
    return "{}/images_{}_{:05d}-of-{:05d}.tfrecord".format(
        output_dir, split_name, shard_id + 1, NUM_SHARDS
    )


def _convert_dataset(split_name, filenames, tps, Qs, classids, output_dir, NUM_SHARDS):
    sys.path.append(os.path.dirname(os.getcwd()))

    from myslim.datasets import dataset_utils

    """Converts the given filenames to a TFRecord dataset.
    Args:
      split_name: The name of the dataset, either 'train' or 'validation'.
      filenames: A list of absolute paths to png or jpg images.
      class_names_to_ids: A dictionary from class names (strings) to ids
                    (integers).
      dataset_dir: The directory where the converted datasets are stored.
    """
    assert split_name in ["train", "validation"]

    num_per_shard = int(math.ceil(len(filenames) / float(NUM_SHARDS)))

    with tf.Graph().as_default():
        image_reader = ImageReader()

        with tf.Session("") as sess:

            for shard_id in range(NUM_SHARDS):
                output_filename = _get_dataset_filename(
                    output_dir, split_name, shard_id, NUM_SHARDS
                )

                print(output_filename)
                with tf.python_io.TFRecordWriter(output_filename) as tfrecord_writer:
                    start_ndx = shard_id * num_per_shard
                    end_ndx = min((shard_id + 1) * num_per_shard, len(filenames))
                    for i in range(start_ndx, end_ndx):
                        # Give status update
                        sys.stdout.write(
                            "\r>> Converting image %d/%d shard %d "
                            % (i + 1, len(filenames), shard_id)
                        )
                        sys.stdout.flush()

                        # Read the filename:
                        print("reading file" + filenames[i])
                        image_data = tf.gfile.FastGFile(filenames[i], "rb").read()
                        height, width = image_reader.read_image_dims(sess, image_data)
                        class_id = classids[i]
                        tp = tps[i]
                        filename = filenames[i]
                        Q = Qs[i]
                        example = dataset_utils.image_to_tfexample(
                            image_data,
                            b"jpg",
                            height,
                            width,
                            int(class_id),
                            int(tp),
                            filename,
                            Q,
                        )
                        tfrecord_writer.write(example.SerializeToString())

            sys.stdout.write("\n")
            sys.stdout.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--slides_folder", help="Set slides folder")
    parser.add_argument("--output_folder", help="Set output folder")
    parser.add_argument("--N_shards", help="Number of shards", default=320)
    args = parser.parse_args()

    output_folder = args.output_folder

    file_info_path = output_folder + "/file_info_train.txt"
    process_train_folder = output_folder + "/process_train"
    if not os.path.exists(process_train_folder):
        os.makedirs(process_train_folder)
    file_info = pd.read_csv(file_info_path, sep="\t")

    training_filenames = list(file_info["tile_path"].values)
    training_classids = [int(id) for id in list(file_info["class_id"].values)]
    tps = [int(id) for id in list(file_info["percent_tumor_cells"].values)]
    Qs = list(file_info["jpeg_quality"].values)

    # First, convert the training and validation sets.
    _convert_dataset(
        split_name="train",
        filenames=training_filenames,
        tps=tps,
        Qs=Qs,
        classidds=training_classids,
        output_dir=process_train_folder,
        NUM_SHARDs=args.N_shards,
    )

    print("Finished converting dataset")
    print(f"The converted data is stored in the directory: {process_train_folder}")
