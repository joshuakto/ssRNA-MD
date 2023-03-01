import subprocess
import os
import numpy as np
import tensorflow.compat.v2 as tf
tf.enable_v2_behavior()
import tensorflow_datasets as tfds
import tensorflow_probability as tfp
tfk = tf.keras
tfkl = tf.keras.layers
tfpl = tfp.layers
tfd = tfp.distributions

def _preprocess(sample):
  ndarray = tf.cast(sample, tf.float32) / 255.  # Scale to unit interval.
  return ndarray, ndarray

def define_model(input_shape, encoded_size, base_depth, plot_model=True):
    prior = tfd.Independent(tfd.Normal(loc=tf.zeros(encoded_size), scale=1),
                        reinterpreted_batch_ndims=1)
    encoder = tfk.Sequential([
        tfkl.InputLayer(input_shape=input_shape),
        tfkl.Lambda(lambda x: tf.cast(x, tf.float32) - 0.5),
        tfkl.Conv2D(base_depth, 5, strides=1,
                    padding='same', activation=tf.nn.leaky_relu),
        tfkl.Conv2D(base_depth, 5, strides=2,
                    padding='same', activation=tf.nn.leaky_relu),
        tfkl.Conv2D(2 * base_depth, 5, strides=1,
                    padding='same', activation=tf.nn.leaky_relu),
        tfkl.Conv2D(2 * base_depth, 5, strides=2,
                    padding='same', activation=tf.nn.leaky_relu),
        tfkl.Conv2D(4 * encoded_size, 7, strides=1,
                    padding='valid', activation=tf.nn.leaky_relu),
        tfkl.Flatten(),
        tfkl.Dense(tfpl.MultivariateNormalTriL.params_size(encoded_size),
                   activation=None),
        tfpl.MultivariateNormalTriL(
            encoded_size,
            activity_regularizer=tfpl.KLDivergenceRegularizer(prior)),
    ])
    decoder = tfk.Sequential([
        tfkl.InputLayer(input_shape=[encoded_size]),
        tfkl.Reshape([1, 1, encoded_size]),
        tfkl.Conv2DTranspose(2 * base_depth, 25, strides=1,
                             padding='valid', activation=tf.nn.leaky_relu),
        tfkl.Conv2DTranspose(2 * base_depth, 5, strides=1,
                             padding='same', activation=tf.nn.leaky_relu),
        tfkl.Conv2DTranspose(2 * base_depth, 5, strides=2,
                             padding='same', activation=tf.nn.leaky_relu),
        tfkl.Conv2DTranspose(base_depth, 5, strides=1,
                             padding='same', activation=tf.nn.leaky_relu),
        tfkl.Conv2DTranspose(base_depth, 5, strides=2,
                             padding='same', activation=tf.nn.leaky_relu),
        tfkl.Conv2DTranspose(base_depth, 5, strides=1,
                             padding='same', activation=tf.nn.leaky_relu),
        tfkl.Conv2D(filters=1, kernel_size=5, strides=1,
                    padding='same', activation=None),
        tfkl.Flatten(),
        tfpl.IndependentBernoulli(input_shape, tfd.Bernoulli.logits),
    ])
    vae = tfk.Model(inputs=encoder.inputs,
                    outputs=decoder(encoder.outputs[0]))

    if plot_model:
        tf.keras.utils.plot_model(vae, show_shapes=True)
        tf.keras.utils.plot_model(decoder, show_shapes=True)
        tf.keras.utils.plot_model(encoder, show_shapes=True)
    return vae, encoder, decoder

def main(dataset_name, encoded_size = 64):
    tf.config.experimental.set_visible_devices([], 'GPU')

    # verify that we have access to a GPU.
    if tf.test.gpu_device_name() != '/device:GPU:0':
      print('WARNING: GPU device not found.')
    else:
      print('SUCCESS: Found GPU: {}'.format(tf.test.gpu_device_name()))

    dataset_name = "distributed100atoms_dataset"
    train_dataset = tf.data.Dataset.load("../dataset/m1_ssRNA_train_" + dataset_name)
    test_dataset = tf.data.Dataset.load("../dataset/m1_ssRNA_test_" + dataset_name)
    input_shape = train_dataset.element_spec[0].shape
    train_dataset = (train_dataset
                     .shuffle(int(10e3))
                     .map(lambda data, label: data)
                     .map(_preprocess)
                     .batch(5)
                     .prefetch(tf.data.AUTOTUNE))
    test_dataset = (test_dataset
                    .shuffle(int(10e3))
                    .map(lambda data, label: data)
                    .map(_preprocess)
                    .batch(5)
                    .prefetch(tf.data.AUTOTUNE))

    # check size of dataset
    print(f'train: {train_dataset.cardinality().numpy()}')
    print(f'test: {test_dataset.cardinality().numpy()}')
    
    base_depth = 32
    vae, encoder, decoder = define_model(input_shape, encoded_size, base_depth)
    callbacks = [
        tfk.callbacks.TensorBoard(
            log_dir="../Tensorboard_log/m1_" + dataset_name,
            histogram_freq=1,  # How often to log histogram visualizations
            embeddings_freq=1,  # How often to log embedding visualizations
        )  
    ]
    negloglik = lambda x, rv_x: -rv_x.log_prob(x)

    vae.compile(optimizer=tf.optimizers.Adam(learning_rate=1e-3),
                loss=negloglik)

    _ = vae.fit(train_dataset,
                epochs=15,
                callbacks=callbacks,
                validation_data=test_dataset)

    vae.save_weights("trained_model/encoded64_m1_ssRNA" + dataset_name)
    return
    

if __name__ == "__main__":
    dataset_name = "distributed100atoms_dataset"
    main(dataset_name)
    print(f"Training using: {dataset_name} done.")