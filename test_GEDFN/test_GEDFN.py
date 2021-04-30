## use desc_tfv1 environment:
# python=3.6.11
# tensorflow=1.15.0
# keras=2.1.0
import os, sys, math
import anndata
import numpy as np
import pandas as pd
import tensorflow as tf

from sklearn.utils import shuffle
from sklearn.preprocessing import OneHotEncoder
from sklearn import metrics

from random import seed
from numpy.random import seed
from tensorflow import set_random_seed

# GPU settings and reproducible
#os.environ["CUDA_VISIBLE_DEVICES"]="1"

#Set seeds
#RANDOM_SEED=0
#seed(RANDOM_SEED)
#np.random.seed(RANDOM_SEED)
#set_random_seed(RANDOM_SEED)

## ----- GEDFN pipeline
class GEDFN(object):
    '''GEDFN pipeline
    '''
    def __init__(self, dims, partition=None, GE=False, gamma_c=50):
        '''
        @dims: network structure
        @partition: gene network adjacency
        @GE: whether use partition or not
        @gamma_c: constant limit for feature selection
        '''
        self.dims = dims
        self.n_layers = len(dims)
        self.n_celltypes = 0
        self.partition = partition
        self.GE = GE
        self.gamma_c = gamma_c

    def init_GEDFN_MLP_params(self):
        weights, biases = {}, {}
        ## set weights and biases
        ## weights from truncated normal and biases are 0
        weights['h1'] =  tf.Variable(tf.random.truncated_normal(shape=[self.n_features, self.n_features], seed=self.random_state, stddev=0.1))
        biases['b1'] = tf.Variable(tf.zeros([self.n_features]))
        for i in range(self.n_layers):
            if i == 0:
                weights['h'+str(i+2)] = tf.Variable(tf.random.truncated_normal(shape=[self.n_features, self.dims[i]], seed=self.random_state, stddev=0.1))
            else:
                weights['h'+str(i+2)] = tf.Variable(tf.random.truncated_normal(shape=[self.dims[i-1], self.dims[i]], seed=self.random_state, stddev=0.1))
            biases['b'+str(i+2)] = tf.Variable(tf.zeros([self.dims[i]]))
        weights['out'] = tf.Variable(tf.random.truncated_normal(shape=[self.dims[-1], self.n_celltypes], seed=self.random_state, stddev=0.1))
        biases['out'] = tf.Variable(tf.zeros([self.n_celltypes]))

        self.weights = weights
        self.biases = biases

    def GEDFN_MLP(self, x, dropout_rate=0.1):
        ## if embedded GE
        if self.GE == True:
            layer_1 = tf.add(tf.matmul(x, tf.multiply(self.weights['h1'], self.partition)), self.biases['b1'])
        else:
            layer_1 = tf.add(tf.matmul(x, self.weights['h1']), self.biases['b1'])
        layer_1 = tf.nn.relu(layer_1)
        if self.max_pooling:
            layer_1 = max_pool(layer_1)
        if self.droph1:
            layer_1 = tf.nn.dropout(layer_1, rate=dropout_rate)
        
        ## add following MLP layers
        layer_prev = layer_1
        for i in range(self.n_layers):
            layer_cur = tf.add(tf.matmul(layer_prev, self.weights['h'+str(i+2)]), self.biases['b'+str(i+2)])
            layer_cur = tf.nn.relu(layer_cur)
            layer_cur = tf.nn.dropout(layer_cur, rate=dropout_rate, seed=self.random_state)
            layer_prev = layer_cur
        out_layer = tf.matmul(layer_prev, self.weights['out']) + self.biases['out']
        return out_layer

    def train(self, x_train, y_train, tf_session=None, batch_size=16, learning_rate=0.0001,
            max_epochs=1000, display_step=1, L2=False, max_pooling=False,
            droph1=False, dropout_rate=0.1, optimizer="Adam", random_state=0,
            avg_cost_criteria=0.1):
        '''
        @x_train: n*p
        @y_train: n*n_celltypes
        '''
        ## initialize network parameters
        self.max_pooling = max_pooling
        self.droph1 = droph1

        ## get number of cell types from y_train
        self.n_celltypes = len(set(y_train.argmax(1)))
        if self.n_celltypes == 0:
            raise ValueError("There should be at least one cell type.")
        self.n_features = x_train.shape[1]
        self.random_state = random_state

        ## tf placeholder for train
        self.x_holder = tf.compat.v1.placeholder(tf.float32, [None, self.n_features], name="InputData")
        self.y_holder = tf.compat.v1.placeholder(tf.int32, [None, self.n_celltypes], name="LabelData")
        self.rate_holder = tf.compat.v1.placeholder(tf.float32)

        ## construct model
        self.init_GEDFN_MLP_params()
        self.pred = self.GEDFN_MLP(self.x_holder, dropout_rate=dropout_rate)

        ## Define loss and optimizer
        tf.stop_gradient(self.y_holder)  # stop gradients flowing on labels
        cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=self.pred, labels=self.y_holder))
        if L2:
            reg = 0
            for key in weights:
                reg += tf.nn.l2_loss(weights[key])
            cost = tf.reduce_mean(cost + 0.01 * reg)

        if optimizer == "Adam":
            self.optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)

        ## Evaluation
        correct_pred = tf.equal(tf.argmax(self.pred, 1), tf.argmax(self.y_holder, 1))
        accuracy = tf.reduce_mean(tf.cast(correct_pred, tf.float32), name="accuracy")
        y_score = tf.nn.softmax(logits=self.pred)

        ## initiate training logs
        loss_rec = np.zeros([max_epochs, 1])
        training_eval = np.zeros([max_epochs, 3])  # include accuracy, ARI and macroF1

        ## open a tf session
        sess = tf_session
        if tf_session is None:
            raise ValueError("Please make sure to initiate a Tensorflow Session before training!")

        ## init variables
        init = tf.compat.v1.global_variables_initializer()
        sess.run(init)
        total_batch = int(np.shape(x_train)[0] / batch_size)

        ## Training cycle
        for epoch in range(max_epochs):
            avg_cost = 0.
            x_tmp, y_tmp = shuffle(x_train, y_train, random_state=random_state)  ## for reproducibility
            # Loop over all batches
            for i in range(total_batch-1):
                batch_x, batch_y = x_tmp[i*batch_size:i*batch_size+batch_size], \
                                    y_tmp[i*batch_size:i*batch_size+batch_size]

                # use cost to track when to stop
                _, c= sess.run([optimizer, cost], feed_dict={self.x_holder: batch_x, 
                                                             self.y_holder: batch_y,
                                                             self.rate_holder: dropout_rate
                                                             })
                # Compute average loss
                avg_cost += c / total_batch

            del x_tmp
            del y_tmp

            ## Display logs per epoch step
            if epoch % display_step == 0:
                loss_rec[epoch] = avg_cost
                acc, y_s = sess.run([accuracy, y_score], 
                        feed_dict={self.x_holder: x_train, 
                                   self.y_holder: y_train, 
                                   self.rate_holder: 0
                                   })

                ## evaluation metrics
                ARI = metrics.cluster.adjusted_rand_score(y_train.argmax(1), y_s.argmax(1))
                macroF1 = metrics.f1_score(y_train.argmax(1), y_s.argmax(1), average="macro")
                training_eval[epoch] = [acc, ARI, macroF1]
                print ("Epoch:", '%d' % (epoch+1), "cost =", "{:.9f}".format(avg_cost),
                        "Training accuracy:", round(acc, 3), " Training ARI:", round(ARI ,3),
                        "Training MacroF1:", round(macroF1, 3))

            if avg_cost <= avg_cost_criteria:
                print("Early stopping when average cost is less than " + str(avg_cost_criteria))

                ## delete unused entries in loss_rec and training_eval
                loss_rec = loss_rec[range(epoch+1)]
                training_eval = training_eval[range(epoch+1)]
                break

        self.loss_rec = loss_rec
        self.training_eval = training_eval

    def predict(self, x_test, tf_session):
        y_score = tf.nn.softmax(logits=self.pred)
        ## Testing cycle
        y_s = tf_session.run([y_score], 
                feed_dict={self.x_holder: x_test, 
                           self.rate_holder: 0
                           }) 
        return y_s[0]

    def feature_importance(self, tf_session):
        if self.GE == True and self.partition is not None:
            ## get features from partition
            gamma_numerator = np.sum(self.partition, axis=0)
            gamma_denominator = np.sum(self.partition, axis=0)
            gamma_numerator[np.where(gamma_numerator>self.gamma_c)] = self.gamma_c
            ## find important features
            var_left = tf.reduce_sum(tf.abs(tf.multiply(self.weights['h1'], self.partition)), 0)
            var_right = tf.reduce_sum(tf.abs(self.weights['h2']), 1)
            var_importance = tf.add(tf.multiply(tf.multiply(var_left, gamma_numerator), 1./gamma_denominator), var_right)

            var_imp = tf_session.run([var_importance])
            var_imp = np.reshape(var_imp, [self.n_features])
            return var_imp
        return None
 
