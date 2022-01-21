#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 16:57:38 2019
@author: arjun
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 12:11:18 2019

@author: Arjun
"""

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Conv2D, MaxPooling2D, Activation, Flatten
from tensorflow.keras import optimizers
from tensorflow.keras.callbacks import TensorBoard, EarlyStopping, ModelCheckpoint
import pickle
import time
import os

print('getcwd:      ', os.getcwd())

X = pickle.load(open("X_data_20sensor_1layer.pickle", "rb"))  # load feature data
y = pickle.load(open("Y_data_20sensor_1layer.pickle", "rb"))  # load target data

print(X.shape)
print(y[1:10])
dense_layers = [0, 1, 3, 4]
layer_sizes_dense = [64, 128, 256]
layer_sizes_convolution = [128, 256, 378, 512]
convolution_layers_1 = [1, 2, 3]
convolution_layers_2 = [1, 2, 3]
convolution_layers_3 = [1, 2, 3]
EPOCHS = 100
batch_size = 256

for dense_layer in dense_layers:
    for size_layer_dense in layer_sizes_dense:
        for convolution_layer_1 in convolution_layers_1:
            for convolution_layer_2 in convolution_layers_2:
                for convolution_layer_3 in convolution_layers_3:
                    for size_layer_convolution in layer_sizes_convolution:

                        model_architecture = "{}denselayer-{}layer_size_dense-{}convlayer1-{}convlayer2-{}convlayer3-{}layer_size_conv-{}".format(
                            dense_layer, size_layer_dense, convolution_layer_1, convolution_layer_2,
                            convolution_layer_3, size_layer_convolution, int(time.time()))
                        print(model_architecture)
                        tensorboard = TensorBoard(log_dir='logs/{}'.format(model_architecture))
                        checkpoint_path = model_architecture
                        checkpoint_dir = os.path.dirname(checkpoint_path)

                        EarlyStopping(monitor='val_loss', patience=10, verbose=1)
                        ModelCheckpoint(filepath=checkpoint_path,
                                           save_weights_only=True,
                                           verbose=1)

                        model = Sequential()

                        # First Group of Conv Layers
                        model.add(Conv2D(size_layer_convolution, (3, 3), input_shape=X.shape[1:], padding='valid'))
                        model.add(Activation("relu"))
                        # model.add(MaxPooling2D(pool_size=(2,2)))

                        for l in range(convolution_layer_1 - 1):
                            model.add(Conv2D(size_layer_convolution, (3, 3), padding='valid'))
                            model.add(Activation("relu"))

                        model.add(MaxPooling2D(pool_size=(2, 2)))

                        # Second Group of Conv Layers
                        for k in range(convolution_layer_2):
                            model.add(Conv2D(size_layer_convolution, (3, 3), padding='valid'))
                            model.add(Activation("relu"))

                        model.add(MaxPooling2D(pool_size=(2, 2)))

                        # Third Group of Conv Layers
                        for e in range(convolution_layer_3):
                            model.add(Conv2D(size_layer_convolution, (3, 3), padding='valid'))
                            model.add(Activation("relu"))

                        model.add(MaxPooling2D(pool_size=(2, 2)))

                        # Fully Connected Layers

                        for _ in range(dense_layer):
                            model.add(Dense(size_layer_dense))
                            model.add(Activation("relu"))
                            model.add(Dropout(0.2))

                        model.add(Flatten())
                        model.add(Dense(1))
                        # model.add(Activation("linear"))
                        model.summary()

                        adam_opt = optimizers.Adam(lr=0.0001, decay=1e-6)
                        model.compile(loss="mean_absolute_error", optimizer=adam_opt, metrics=['mae', 'mse'])
                        model.fit(X, y, batch_size=batch_size, epochs=EPOCHS, validation_split=0.3, callbacks=[tensorboard])

