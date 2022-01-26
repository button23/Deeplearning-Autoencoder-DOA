import datetime
import random
import time
import numpy as np
import tensorflow as tf
import random as rn
import os
import json

from tensorflow.keras.optimizers import Adam
# from keras.optimizers import metrics, regularizers, optimizers, backend optimizers

from tensorflow.keras.callbacks import TensorBoard, EarlyStopping, ModelCheckpoint
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization, Activation, Conv1D, Flatten, MaxPool1D, Conv2D, ZeroPadding1D, concatenate
from tensorflow.keras.utils import to_categorical
import scipy.io as sio
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score

dense_layers = [0, 1, 2, 3]
dense_sizes = [32, 64, 128, 256, 512]
conv_layers = [1, 2, 3, 4]
num_filters = [32, 64, 128, 256]
kernel_sizes = [2000, 2250, 2500, 2750, 3000]

# number of reference point
nClass = 180
nSample = 56
# Model configuration
num_epochs = 40
batch_size = 128
validation_split = 0.2
rootPath = os.getcwd() + '/data/1224/train/'

# gpu setting
# gpu_options = tf.compat.v1.GPUOptions(per_process_gpu_memory_fraction=0.333)
# sess = tf.compat.v1.Session(config=tf.compat.v1.ConfigProto(gpu_options=gpu_options))

dataPath = rootPath + 'train_data.mat'
labelPath = rootPath + 'train_label.mat'
x_train = sio.loadmat(dataPath)
y_train = sio.loadmat(labelPath)

rootPath = os.getcwd() + '/data/1224/test/'
dataPath = rootPath + 'test_data.mat'
labelPath = rootPath + 'test_label.mat'
x_test = sio.loadmat(dataPath)
y_test = sio.loadmat(labelPath)

# load data
x_train = x_train['train_data']
y_train = y_train['train_label']
print(x_train)

x_test = x_test['test_data']
y_test = y_test['test_label']
print(x_test)

x_train_real = x_train.real
x_train_imag = x_train.imag
x_train = np.concatenate((x_train_real, x_train_imag), axis=-1)
print(x_train.shape)

# create scaler (Normalization)
scaler = MinMaxScaler()
# scaler = StandardScaler()
x_train = scaler.fit_transform(x_train)
# # hot encoding
y_train = to_categorical(y_train, nClass)
# train data shuffle
np.random.seed(2020)
index = np.arange(x_train.shape[0])
print(index)
np.random.shuffle(index)
x_train = x_train[index]
y_train = y_train[index]
x_train = x_train.reshape(x_train.shape[0], x_train.shape[1], 1)


ker_size = 3
# for dense_layer in dense_layers:
#     for conv_layer in conv_layers:
#         ker_size = random.choice(kernel_sizes)
#         num_filter = random.choice(num_filters)
        # define multiple inputs
_input_ = Input(shape=(x_train.shape[1], 1))

_hid_ = Conv1D(filters=16, kernel_size=ker_size, strides=1, padding='valid', use_bias=True,
               kernel_initializer='he_uniform')(_input_)
_hid_ = BatchNormalization(epsilon=1e-06, momentum=0.9, weights=None)(_hid_)
_hid_ = Activation('relu')(_hid_)
_hid_ = MaxPool1D(pool_size=ker_size, strides=1, padding='same', data_format=None)(_hid_)

# for l in range(conv_layer - 1):
#     if l != 0:
_hid_ = Conv1D(filters=32, kernel_size=ker_size, strides=1, padding='valid', use_bias=True,
               kernel_initializer='he_uniform')(_hid_)
_hid_ = BatchNormalization(epsilon=1e-06, momentum=0.9, weights=None)(_hid_)
_hid_ = Activation('relu')(_hid_)
_hid_ = MaxPool1D(pool_size=ker_size, strides=1, padding='same', data_format=None)(_hid_)

_hid_ = Flatten()(_hid_)
# for l in range(dense_layer):
#     dense_size = random.choice(dense_sizes)
_hid_ = Dense(32, activation='relu', kernel_initializer='he_uniform')(_hid_)

_output_ = Dense(nClass, activation='softmax', kernel_initializer='he_uniform')(_hid_)

model = Model(inputs=_input_, outputs=_output_)
# model.summary()

# NAME = "{}-conv-{}-kersize-{}-dense-{}".format(conv_layer, ker_size, dense_layer, int(time.time()))
NAME = "2-conv-16-32-kersize-{}-1-dense-32-{}".format(ker_size, int(time.time()))

savePath =rootPath + 'mobility_model/{}.h5'.format(NAME)  # 256_1024conv_128_64dense.h5

# callbacks_list
callbacks_list = [
    EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=30),
    # EarlyStopping(monitor='val_loss', patience=20),
    ModelCheckpoint(filepath=savePath, monitor='val_loss', verbose=0, save_best_only=True,
                    save_weights_only=False,
                    mode='min', period=1),
    TensorBoard(log_dir='./logs/{}'.format(NAME), histogram_freq=1, write_graph=False, write_images=False)
]

# Optimizser
adam = Adam(lr=1e-3, beta_1=0.9, beta_2=0.999, epsilon=1e-08, amsgrad=True, decay=1e-5)

# compile and fit
model.compile(loss='categorical_crossentropy', optimizer=adam, metrics=['categorical_accuracy'])
# model.fit([x_train0, x_train1, x_train2], y_train, epochs=num_epochs, batch_size=batch_size,
model.fit(x_train, y_train, epochs=num_epochs, batch_size=batch_size,
          validation_split=validation_split, shuffle=True, callbacks=callbacks_list)

            # evaluate
            # scores = model.evaluate(x_train, y_train)
            # print(" %s %f" % (model.metrics_names[1], scores[1]))


def predict(model, data):
  reconstructions = model(data)
  return reconstructions

def print_stats(predictions, labels):
  print("Accuracy = {}".format(accuracy_score(labels, predictions)))
  print("Precision = {}".format(precision_score(labels, predictions)))
  print("Recall = {}".format(recall_score(labels, predictions)))

preds = predict(model, x_test)
print_stats(preds, test_labels)