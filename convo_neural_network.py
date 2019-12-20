import tensorflow as tf
import numpy as np
from tensorflow.keras import datasets, layers, datasets, models
features = np.load('./Data/Clean/features.npy')
features = np.expand_dims(features, 2)
print(features.shape)
outputs = np.load('./Data/Clean/outputs.npy')
n_data = features.shape[0]
n_features = features.shape[1]
model = models.Sequential()
model.add(layers.Conv1D(filters=1, kernel_size=100, activation='tanh', input_shape=(n_features,1)))
model.add(layers.MaxPooling1D())
model.add(layers.Flatten())
model.add(layers.Dense(units=1000, activation='tanh'))
model.add(layers.Dense(units=100, activation='tanh'))
model.add(layers.Dense(units=1))
model.compile(loss='mean_squared_error', optimizer='adam', metriocs=['accuracy'])
model.fit(features, outputs, epochs=5, batch_size=50)



