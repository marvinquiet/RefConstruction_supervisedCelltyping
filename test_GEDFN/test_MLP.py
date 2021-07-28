# TensorFlow and tf.keras
import tensorflow.compat.v1 as tf
from tensorflow.python.ops import array_ops, init_ops

from tensorflow import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout

# Helper libraries
import numpy as np
import matplotlib.pyplot as plt

#print(tf.__version__)

class MLP(object):
    def __init__(self, dims):
        self.dims = dims

    def focal_loss(self, y_true, y_pred, gamma=0, alpha=1):
        '''Adapted from tensorflow addons, but also refer to other links
        #https://github.com/fizyr/keras-retinanet/blob/424671f71da40845e987713544493f4e88f4ea68/keras_retinanet/losses.py
        #https://github.com/ailias/Focal-Loss-implement-on-Tensorflow/blob/master/focal_loss.py
        #https://github.com/umbertogriffo/focal-loss-keras/blob/master/src/loss_function/losses.py
        https://github.com/tensorflow/addons/blob/11aad761bfdc5cd4ef29b29e678e0d081f846fa3/tensorflow_addons/losses/focal_loss.py
        https://github.com/artemmavrin/focal-loss/blob/5b2ca68/src/focal_loss/_categorical_focal_loss.py

        .. math::
            L(y, \hat{\mathbf{p}})
                = -\left(1 - \hat{p}_y\right)^\gamma \log(\hat{p}_y)

        Focal loss for imbalanced classification
        
        @y_true: target data
        @y_pred: predicted data
        @gamma: modulating factor, change the focus from well-classified to hard classes
        @alpha: balancing factor, either 1 or None (inverse frequency, or inverse document frequency)
        '''
        if None == alpha:
            alpha = self.weights

        alpha = np.array(alpha, dtype=np.float32)

        epsilon = keras.backend.epsilon()
        y_pred = keras.backend.clip(y_pred, epsilon, 1.-epsilon)
        ce = -y_true * keras.backend.log(y_pred)

        loss = alpha * keras.backend.pow(1-y_pred, gamma) * ce
        return keras.backend.mean(keras.backend.sum(loss, axis=-1))

    def init_MLP(self):
        dense_kernel_init = keras.initializers.TruncatedNormal(mean=0, stddev=0.1, seed=self.random_state) ## same as GEDFN
        model = Sequential()
        for i in range(len(self.dims)):
            if 0 == i:
                model.add(Dense(self.dims[i], input_shape=self.input_shape, 
                    activation=tf.nn.relu,
                    kernel_initializer=dense_kernel_init))
            else:
                model.add(Dense(self.dims[i], activation=tf.nn.relu,
                    kernel_initializer=dense_kernel_init))
            model.add(Dropout(rate=self.dropout_rate, input_shape=(self.dims[i],), seed=self.random_state))
        model.add(Dense(self.n_classes, activation=tf.nn.softmax,
            kernel_initializer=dense_kernel_init))
        #model.compile(loss=tf.keras.losses.CategoricalCrossentropy(), 
        model.compile(loss=[self.focal_loss], 
                metrics=["accuracy"], ## show training accuracy
                optimizer=self.optimizer)
        return model
    
    def fit(self, x_train, y_train, batch_size=16, max_epochs=1000, dropout_rate=0.1, 
            optimizer="adam", random_state=0):
        ## get input shape and initialize parameters
        self.input_shape = (x_train.shape[1], )
        self.n_classes = len(set(y_train.argmax(1)))

        ## calculate weights for each sample using inverse frequency
        self.weights = inverse_variance_weighting(y_train.argmax(1))

        print(self.weights)

        self.dropout_rate = dropout_rate
        self.optimizer = optimizer
        self.random_state = random_state

        ## init MLP model
        self.model = self.init_MLP()
        #print(self.model.summary())

        ## add callback with 5 steps no improvement
        callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=5)
        #self.model.fit(x_train, y_train, epochs=max_epochs, batch_size=batch_size, 
        #        validation_split=0.2, callbacks=[callback], verbose=2)

        self.model.fit(x_train, y_train, epochs=max_epochs, batch_size=batch_size, 
                validation_split=0.0, callbacks=[callback], verbose=2) # without cross validation

    def predict(self, x_test):
        x_pred = self.model.predict(x_test)
        return x_pred

def inverse_frequency(logits):
    ''' The inverse frequency of labels

    @logits: a list of labels
    '''
    unique, counts = np.unique(logits, return_counts=True)
    print(unique, counts)
    freq = counts/sum(counts)  ## get frequency
    inv_freq = 1/freq   ## get inverse frequency

    return inv_freq/sum(inv_freq) ## normalize to 1


def inverse_document_frequency(logits):
    '''Use inverse document frequency with smooth
    '''
    unique, counts = np.unique(logits, return_counts=True)
    print(unique, counts)

    freq = np.log(sum(counts)/(1+counts)) + 1
    return freq/sum(freq)  ## normalize to 1

def inverse_variance_weighting(logits):
    '''Inverse variance weighting, using square root
    '''
    unique, counts = np.unique(logits, return_counts=True)
    print(unique, counts)
    freq = counts/sum(counts)
    inv_freq = np.sqrt(1/freq)  ## square root inverse frequency

    return inv_freq/sum(inv_freq) ## normalize to 1
