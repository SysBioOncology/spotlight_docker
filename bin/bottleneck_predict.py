#!/usr/bin/env python
import os
import sys
import time

import tensorflow.compat.v1 as tf
import tf_slim as slim
from nets import nets_factory
from preprocessing import preprocessing_factory

tf.compat.v1.disable_eager_execution()

tf.app.flags.DEFINE_integer("num_classes", 42, "The number of classes.")
tf.app.flags.DEFINE_string("bot_out", None, "Output file for bottleneck features.")
tf.app.flags.DEFINE_string("pred_out", None, "Output file for predictions.")
tf.app.flags.DEFINE_string(
    "model_name", "inception_v4", "The name of the architecture to evaluate."
)
tf.app.flags.DEFINE_string(
    "checkpoint_path", None, "The directory where the model was written to."
)
tf.app.flags.DEFINE_integer("eval_image_size", 299, "Eval image size.")
tf.app.flags.DEFINE_string("file_dir", "../Output/process_train/", "")

FLAGS = tf.app.flags.FLAGS


def main(_):
    model_name_to_variables = {
        "inception_v3": "InceptionV3",
        "inception_v4": "InceptionV4",
    }
    model_name_to_bottleneck_tensor_name = {
        "inception_v4": "InceptionV4/Logits/AvgPool_1a/AvgPool:0",
        "inception_v3": "InceptionV3/Logits/AvgPool_1a_8x8/AvgPool:0",
    }
    bottleneck_tensor_name = model_name_to_bottleneck_tensor_name.get(FLAGS.model_name)
    preprocessing_name = FLAGS.model_name
    eval_image_size = FLAGS.eval_image_size
    model_variables = model_name_to_variables.get(FLAGS.model_name)
    if model_variables is None:
        tf.logging.error("Unknown model_name provided `%s`." % FLAGS.model_name)
        sys.exit(-1)
    # Either specify a checkpoint_path directly or find the path
    if tf.gfile.IsDirectory(FLAGS.checkpoint_path):
        checkpoint_path = tf.train.latest_checkpoint(FLAGS.checkpoint_path)
        print(checkpoint_path)
        if checkpoint_path is None:
            sys.exit(-1)
    else:
        checkpoint_path = FLAGS.checkpoint_path

    image_string = tf.placeholder(tf.string)
    image = tf.image.decode_jpeg(
        image_string, channels=3, try_recover_truncated=True, acceptable_fraction=0.3
    )
    image_preprocessing_fn = preprocessing_factory.get_preprocessing(
        preprocessing_name, is_training=False
    )
    network_fn = nets_factory.get_network_fn(
        FLAGS.model_name, FLAGS.num_classes, is_training=False
    )
    processed_image = image_preprocessing_fn(image, eval_image_size, eval_image_size)
    processed_images = tf.expand_dims(processed_image, 0)

    logits, _ = network_fn(processed_images)
    probabilities = tf.nn.softmax(logits)
    init_fn = slim.assign_from_checkpoint_fn(
        checkpoint_path, slim.get_model_variables(model_variables)
    )

    print(FLAGS.bot_out)

    sess = tf.Session()
    init_fn(sess)

    fto_bot = open(FLAGS.bot_out, "w")
    fto_pred = open(FLAGS.pred_out, "w")

    filelist = [
        file_path
        for file_path in os.listdir(FLAGS.file_dir)
        if (file_path.startswith("images_train") & file_path.endswith(".tfrecord"))
    ]
    for i in range(len(filelist)):
        file = filelist[i]
        fls = tf.python_io.tf_record_iterator(FLAGS.file_dir + "/" + file)
        tf.logging.info("reading from: %s" % file)
        start_time = time.time()
        c = 0
        for fl in fls:
            example = tf.train.Example()
            example.ParseFromString(fl)
            x = example.features.feature["image/encoded"].bytes_list.value[0]
            filenames = str(
                example.features.feature["image/filename"].bytes_list.value[0]
            )
            label = str(
                example.features.feature["image/class/label"].int64_list.value[0]
            )
            preds = sess.run(probabilities, feed_dict={image_string: x})
            bottleneck_values = sess.run(bottleneck_tensor_name, {image_string: x})
            fto_pred.write(filenames + "\t" + label)
            fto_bot.write(filenames + "\t" + label)
            for p in range(len(preds[0])):
                fto_pred.write("\t" + str(preds[0][p]))
            fto_pred.write("\n")
            for p in range(len(bottleneck_values[0][0][0])):
                fto_bot.write("\t" + str(bottleneck_values[0][0][0][p]))
            fto_bot.write("\n")
            c += 1
        used_time = time.time() - start_time
        tf.logging.info("processed images: %s" % c)
        tf.logging.info("used time: %s" % used_time)

    fto_bot.close()
    fto_pred.close()
    sess.close()


if __name__ == "__main__":
    tf.app.run()
