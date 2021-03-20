import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tensorflow as tf
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.initializers import RandomNormal
from tensorflow.keras.layers import Activation
from tensorflow.keras.callbacks import EarlyStopping

# import tensorflow_addons as tfa


def build_model(vocab_size, num_units, dropout, optimizer):

    es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=50)

    weight_init = RandomNormal(mean=0.0, stddev=0.05, seed=71)

    model = tf.keras.Sequential([
        tf.keras.layers.Embedding(vocab_size, num_units), #, input_length=max_length),
        tf.keras.layers.LSTM(num_units, return_sequences=True, kernel_initializer=weight_init, dropout=dropout),
        tf.keras.layers.LSTM(num_units, kernel_initializer=weight_init, dropout=dropout),
        tf.keras.layers.Dense(vocab_size, activation="softmax")
    ])

    model.compile(loss='categorical_crossentropy', optimizer=optimizer, metrics=['mae', 'acc'])

    return model
