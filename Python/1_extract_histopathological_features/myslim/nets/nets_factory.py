# Copyright 2016 The TensorFlow Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================
"""Contains a factory for building various models."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import functools

# Originally
# import tensorflow as tf 
# slim = tf.contrib.slim

# Replaced with
import tensorflow.compat.v1 as tf
import tf_slim as slim

from nets import inception
from nets import inception_v4_alt
from nets import overfeat

networks_map = {'overfeat': overfeat.overfeat,
                'inception_v4': inception.inception_v4,
                #                'inception_v4_alt': inception.inception_v4_alt,
                }

arg_scopes_map = {'overfeat': overfeat.overfeat_arg_scope,
                  'inception_v4': inception.inception_v4_arg_scope,
                  #                  'inception_v4_alt': inception.inception_v4_alt_arg_scope,
                  }


def get_network_fn(name, num_classes, weight_decay=0.0, is_training=False):
    """Returns a network_fn such as `logits, end_points = network_fn(images)`.

    Args:
      name: The name of the network.
      num_classes: The number of classes to use for classification.
      weight_decay: The l2 coefficient for the model weights.
      is_training: `True` if the model is being used for training and `False`
        otherwise.

    Returns:
      network_fn: A function that applies the model to a batch of images. It has
        the following signature:
          logits, end_points = network_fn(images)
    Raises:
      ValueError: If network `name` is not recognized.
    """
    if name not in networks_map:
        raise ValueError('Name of network unknown %s' % name)
    arg_scope = arg_scopes_map[name](weight_decay=weight_decay)
    func = networks_map[name]

    @functools.wraps(func)
    def network_fn(images):
        with slim.arg_scope(arg_scope):
            return func(images, num_classes, is_training=is_training)
    if hasattr(func, 'default_image_size'):
        network_fn.default_image_size = func.default_image_size

    return network_fn
