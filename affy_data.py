#import pyaffy
import pandas as pd
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam

from sklearn.preprocessing import MinMaxScaler
import numpy as np

scaler = MinMaxScaler()


values = pd.read_csv('./Data/Affymetrix/affy_data.csv')
values = values.drop(axis=1, labels=values.columns[0])
values = values.to_numpy()
values = scaler.fit_transform(values)
scaler = MinMaxScaler()
values = values.T

# Autoencoder decoder
input_data = Input(shape=(values.shape[1],))
layer1 = Dense(5000, activation='relu', kernel_initializer='glorot_normal')(input_data)
layer2 = Dense(2000, activation = 'relu', kernel_initializer='glorot_normal')(layer1)
encoded_layer = Dense(1000, activation = 'relu', kernel_initializer='glorot_normal')(layer2)
layer3 = Dense(2000, activation='relu', kernel_initializer='glorot_normal')(encoded_layer)
layer4 = Dense(5000, activation = 'relu', kernel_initializer='glorot_normal')(layer3)
decoded_layer = Dense(values.shape[1], activation='sigmoid', kernel_initializer='glorot_normal')(layer4)

autoencoder = Model(input_data, decoded_layer)
adam_opt = Adam(learning_rate=0.01)
autoencoder.compile(optimizer=adam_opt, loss='binary_crossentropy')
autoencoder.fit(values, values, epochs=85, batch_size=100, shuffle=True)
autoencoder.save('autoencoder_model.h5')