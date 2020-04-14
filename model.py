import tensorflow as tf
import numpy as np
import os
    
class WGAN_GP(object):
    def __init__(self, data, nte):
        #hyperparameters
        self.lr = 1e-5
        self.epochs = 200000
        self.batch = 32
        self.zdim = 100
        self.gensize = 250
        self.discsize = 150
        self.data = data
        self.n = self.data.shape[0]
        self.p = self.data.shape[1]
        self.tr_max = np.max(self.data)
        self.nte = nte
        
        self.fdata = self.data.flatten()
        self.sess = tf.InteractiveSession()
        self.build_model()
        
    """
    def data_split(self, test_ratio=0.1):
        xidx = np.arange(len(self.data))
        self.teidx = np.random.choice(len(self.data), int(len(self.data)*test_ratio), False)
        self.tridx = np.setdiff1d(xidx, self.teidx)
        
        self.ntr = len(self.tridx)
        self.nte = len(self.teidx)
        
        dat_tr, dat_te = self.data[self.tridx], self.data[self.teidx]
        self.n = dat_tr.shape[0]
        self.p = dat_tr.shape[1]
        self.tr_max = np.max(dat_tr)
        return dat_tr, dat_te
    """
    def get_z(self, batch):
        return self.fdata[np.random.randint(0,len(self.fdata),size=(batch,self.zdim))]
    
    def get_batch(self):
        idx = np.random.choice(self.n, self.batch, replace=False)
        return self.data[idx]
    
    def build_model(self):
        self.z = tf.placeholder(tf.float32, shape=(None, self.zdim))
        self.x = tf.placeholder(tf.float32, shape=(None, self.p))
    
        self.G = self.generator(self.z)
        self.D_real = self.discriminator(self.x)
        self.D_fake = self.discriminator(self.G, reuse=True)
        dw1, dw2 = self.discriminator(self.x, types=1, reuse=True)
        
        #gradient penalty
        lambda_ = 10
        eps = tf.random_uniform([self.batch,1], minval=0, maxval=self.tr_max)
        x_inter = eps*self.x + (1-eps)*self.G
        grad = tf.gradients(self.discriminator(x_inter, reuse=True), [x_inter])[0]
        dgrad = tf.sqrt(tf.reduce_sum(tf.square(grad), axis=1))
        gp = lambda_*tf.reduce_mean(tf.square(dgrad-1.))
        
        self.obj_d = tf.reduce_mean(self.D_fake) - tf.reduce_mean(self.D_real) + gp
        self.obj_g = -tf.reduce_mean(self.D_fake)
        
        t_vars = tf.trainable_variables()
        self.d_vars = [var for var in t_vars if var.name.startswith('discriminator')]
        self.g_vars = [var for var in t_vars if var.name.startswith('generator')]
        
        self.saver = tf.train.Saver(max_to_keep=self.epochs)
        
    
    def train(self):
        savepath = os.getcwd()+'/weights/'
        self.savepath = savepath
        if not os.path.exists(savepath):
            os.mkdir(savepath)
        #np.savez(savepath+'idx.npz', tridx=self.tridx, teidx=self.teidx)
            
        opt_g = tf.train.AdamOptimizer(learning_rate=self.lr).minimize(self.obj_g, var_list=self.g_vars)
        opt_d = tf.train.AdamOptimizer(learning_rate=self.lr).minimize(self.obj_d, var_list=self.d_vars)
        
        #self.sess = tf.InteractiveSession()
        init = tf.global_variables_initializer().run()
        
        genx = []
        dloss, gloss = [], []
        tmpz = self.get_z(self.nte)
        self.saver.save(self.sess, savepath+'ep0.ckpt') #save initial weight
        for ep in range(self.epochs):
            sname = savepath+'ep'+str(ep+1)+'.ckpt'
            for it in range(self.n//self.batch):
                self.sess.run(opt_d, {self.x:self.get_batch(), self.z:self.get_z(self.batch)})
                self.sess.run(opt_g, {self.x:self.get_batch(), self.z:self.get_z(self.batch)})
            dl, gl = self.sess.run([self.obj_d, self.obj_g], {self.x:self.get_batch(), self.z:self.get_z(self.batch)})
            dloss.append(dl)
            gloss.append(gl)
        
            if (ep+1)%500==0:
                print ('epochs:', ep+1, '- dloss:', dl, ', gloss:', gl)
                self.saver.save(self.sess, sname)
                tmpx = self.sess.run(self.G, {self.z:tmpz})
                genx.append(tmpx)
        self.genx = genx
        np.savez(savepath+'training_results.npz', dloss=dloss, gloss=gloss, genx=self.genx)
        return dloss, gloss, genx
                
        
    def generator(self, z, reuse=False):
        gen_init = 0.3
        with tf.variable_scope('generator', reuse=reuse):
            gen_dense1 = tf.layers.dense(inputs=z, units=self.gensize, kernel_initializer=tf.random_uniform_initializer(-gen_init,gen_init), activation=tf.nn.leaky_relu, name='gen_dense1')
            gen_dense2 = tf.layers.dense(inputs=gen_dense1, units=self.gensize, kernel_initializer=tf.random_uniform_initializer(-gen_init,gen_init), activation=tf.nn.leaky_relu, name='gen_dense2')
            gen_logit = tf.layers.dense(inputs=gen_dense2, units=self.p, kernel_initializer=tf.random_uniform_initializer(-gen_init,gen_init), activation=None, name='gen_logit')
            return gen_logit
        
    def discriminator(self, x, types=0, reuse=False):
        with tf.variable_scope('discriminator', reuse=reuse):
            disc_dense1 = tf.layers.dense(inputs=x, units=self.discsize, activation=tf.nn.leaky_relu, name='disc_dense1')
            disc_dense2 = tf.layers.dense(inputs=disc_dense1, units=self.discsize, activation=tf.nn.leaky_relu, name='disc_dense2')
            disc_logit = tf.layers.dense(inputs=disc_dense2, units=1, activation=None, name='disc_logit')
            if types==0: #return logit
                return disc_logit
            elif types==1: #for weight calculation
                return disc_dense1, disc_dense2
            
    def load_weight(self, loadpath):
        self.saver.restore(self.sess, loadpath+'.ckpt')
    
    def generate_samples(self, z):
        if type(z)==int: #z is int
            tz = self.get_z(z)
            tgen = self.sess.run(self.G, {self.z:tz})
            return tz, tgen
        else: #z is vector
            tz = z
            tgen = self.sess.run(self.G, {self.z:tz})
            return tgen
        
    def get_weights(self, layer):
        return self.sess.run(layer)
        