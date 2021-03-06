import FVF_model_base
import numpy as np
from numpy import sin, cos, tan

class Model(FVF_model_base.Model):
    '''
    Class to run the FVF model for MAC waves
    defines make_A and make_B

    Equations Include:
        R-momentum
        Theta-momentum
        Phi-momentum
        B-divergence (Gauss' Law)
        Theta-Induction
        Phi-induction
        U-divergence (Mass conservation)
        Displacement Equation
    '''

    def make_A(self):
        Nk = self.Nk
        Nl = self.Nl
        E = self.E
        Pm = self.Pm
        N = self.N
        th = self.th
        Br = self.Br
        Bth = self.Bth
        Bph = self.Bph
        ones = np.ones((Nk,Nl))
        index_top_zero_N= np.max(np.where(N[:,0]==0.))
        index_bottom_stratified = index_top_zero_N + 1

        '''
        Creates the A matrix (M*l*x = A*x)
        m: azimuthal fourier mode to compute
        '''

        ################################
        # Momentum Equation ############
        ################################
        # R-momentum
        self.add_gov_equation('rmom', 'vr')
        self.rmom.add_drP('p', C= -1, k_vals=range(index_bottom_stratified,Nk-1))
        # self.rmom.add_drP('p', C=-1, k_vals=[1])
        # self.rmom.add_drP('p', C=-1, k_vals=[1])
        self.rmom.add_term('ur', -N**2, k_vals=range(index_bottom_stratified,Nk-1))
        self.rmom.add_d2_b0('vr', C= E, k_vals=range(index_bottom_stratified,Nk-1))
        self.rmom.add_d2r_th('vth', C= E, k_vals=range(index_bottom_stratified,Nk-1))
        self.rmom.add_d2r_ph('vph', C= E, k_vals=range(index_bottom_stratified,Nk-1))
        ### set pressure to zero within the non-stratified layer:
        self.rmom.add_term('p', ones, k_vals=range(0,index_bottom_stratified))
        # self.rmom.add_term('p', ones, kdiff=1, k_vals=[0])
        self.rmom.add_term('p', ones, k_vals=[Nk-1])
        self.rmom.add_term('p', -ones, kdiff=-1, k_vals=[Nk-1])

        self.A_rows = self.rmom.rows
        self.A_cols = self.rmom.cols
        self.A_vals = self.rmom.vals
        del self.rmom

        # Theta-Momentum
        self.add_gov_equation('tmom', 'vth')
        self.tmom.add_dthP('p', C= -1, k_vals=range(1,Nk-1))
        self.tmom.add_term('vph', 2.0*cos(th), k_vals=range(1,Nk-1))
        self.tmom.add_d2_bd0('vth', C= E, k_vals=range(1,Nk-1))
        self.tmom.add_d2th_r('vr', C= E, k_vals=range(1,Nk-1))
        self.tmom.add_d2th_ph('vph', C= E, k_vals=range(1,Nk-1))
        self.tmom.add_dr_b0('bth', C= Br*E/Pm, k_vals=range(1,Nk-1))
        # self.tmom.add_term('bth', simple_mag_bc_bdr, k_vals=[0])
        self.tmom.add_term('vth', ones, k_vals=[0])
        # self.tmom.add_term('vth', -ones, kdiff=1, k_vals=[0])
        self.tmom.add_term('vth', ones, k_vals=[Nk-1])
        self.tmom.add_term('vth', -ones, kdiff=-1, k_vals=[Nk-1])
        # self.tmom.add_dth('br', C= -Br*E/Pm, k_vals=range(1,Nk-1))
        self.A_rows += self.tmom.rows
        self.A_cols += self.tmom.cols
        self.A_vals += self.tmom.vals
        del self.tmom

        # Phi-Momentum
        self.add_gov_equation('pmom', 'vph')
        self.pmom.add_dphP('p', C= -1, k_vals=range(1,Nk-1))
        self.pmom.add_term('vth', -2.0*cos(th), k_vals=range(1,Nk-1))
        self.pmom.add_d2_bd0('vph', C= E, k_vals=range(1,Nk-1))
        self.pmom.add_d2ph_r('vr', C= E, k_vals=range(1,Nk-1))
        self.pmom.add_d2ph_th('vth', C= E, k_vals=range(1,Nk-1))
        self.pmom.add_dr_b0('bph', C= Br*E/Pm, k_vals=range(1,Nk-1))
        # self.pmom.add_term('bph', simple_mag_bc_bdr, k_vals=[0])
        self.pmom.add_term('vph', ones, k_vals=[0])
        # self.pmom.add_term('vph', -ones, kdiff=1, k_vals=[0])
        self.pmom.add_term('vph', ones, k_vals=[Nk-1])
        self.pmom.add_term('vph', -ones, kdiff=-1, k_vals=[Nk-1])
        # self.pmom.add_dph('br', C= -Br*E/Pm, k_vals=range(1,Nk-1))
        self.A_rows += self.pmom.rows
        self.A_cols += self.pmom.cols
        self.A_vals += self.pmom.vals
        del self.pmom

        ################################
        # Lorentz Equation ##########
        ################################

        # # B-divergence replaces r-lorentz
        # self.add_gov_equation('bdiv', 'br')
        # self.bdiv.add_dr_bd0('br')
        # self.bdiv.add_dth('bth')
        # self.bdiv.add_dph('bph')
        # self.A_rows += self.bdiv.rows
        # self.A_cols += self.bdiv.cols
        # self.A_vals += self.bdiv.vals
        # del self.bdiv

        # theta-Lorentz
        self.add_gov_equation('thlorentz', 'bth')
        self.thlorentz.add_dr_bd0('vth', C= Br, k_vals=range(1,Nk-1))
        self.thlorentz.add_d2('bth', C= E/Pm, k_vals=range(1,Nk-1))
        # self.thlorentz.add_d2th_r('br', C= E/Pm, k_vals=range(1,Nk-1))
        self.thlorentz.add_d2th_ph('bph', C= E/Pm, k_vals=range(1,Nk-1))
        self.thlorentz.add_term('bth', ones, k_vals=[0])
        self.thlorentz.add_term('bth', ones, k_vals=[Nk-1])
        self.A_rows += self.thlorentz.rows
        self.A_cols += self.thlorentz.cols
        self.A_vals += self.thlorentz.vals
        del self.thlorentz

        # phi-Lorentz
        self.add_gov_equation('phlorentz', 'bph')
        self.phlorentz.add_dr_bd0('vph', C= Br, k_vals=range(1,Nk-1))
        self.phlorentz.add_d2('bph', C= E/Pm, k_vals=range(1,Nk-1))
        # self.phlorentz.add_d2ph_r('br', C= E/Pm, k_vals=range(1,Nk-1))
        self.phlorentz.add_d2ph_th('bth', C= E/Pm, k_vals=range(1,Nk-1))
        self.phlorentz.add_term('bph', ones, k_vals=[0])
        self.phlorentz.add_term('bph', ones, k_vals=[Nk-1])
        self.A_rows += self.phlorentz.rows
        self.A_cols += self.phlorentz.cols
        self.A_vals += self.phlorentz.vals
        del self.phlorentz

        # Divergence (Mass Conservation) #########
        self.add_gov_equation('div', 'p')
        self.div.add_dr_b0('vr', k_vals=range(1,Nk-1))
        self.div.add_dth('vth', k_vals=range(1,Nk-1))
        self.div.add_dph('vph', k_vals=range(1,Nk-1))

        self.div.add_term('vr', ones, k_vals=[0])
        self.div.add_term('vr', ones, k_vals=[Nk-1])

        self.A_rows += self.div.rows
        self.A_cols += self.div.cols
        self.A_vals += self.div.vals
        del self.div

        # Displacement Equation #########
        self.add_gov_equation('ur', 'ur')
        self.ur.add_term('vr', ones)
        self.A_rows += self.ur.rows
        self.A_cols += self.ur.cols
        self.A_vals += self.ur.vals
        del self.ur

        self.A = FVF_model_base.coo_matrix((self.A_vals, (self.A_rows, self.A_cols)),
                                   shape=(self.SizeM, self.SizeM))
        del self.A_vals, self.A_rows, self.A_cols
        return self.A

    def make_B(self):
        '''
        Creates the B matrix (B*l*x = A*x)
        m: azimuthal fourier mode to compute
        '''
        ones = np.ones((self.Nk, self.Nl))
        Nk = self.Nk
        Nl = self.Nl
        self.B_rows = []
        self.B_cols = []
        self.B_vals = []

        self.add_gov_equation('B_vth', 'vth')
        self.B_vth.add_term('vth', ones, k_vals=range(1,Nk-1))
        self.B_rows = self.B_vth.rows
        self.B_cols = self.B_vth.cols
        self.B_vals = self.B_vth.vals
        del self.B_vth

        self.add_gov_equation('B_vph', 'vph')
        self.B_vph.add_term('vph', ones, k_vals=range(1,Nk-1))
        self.B_rows += self.B_vph.rows
        self.B_cols += self.B_vph.cols
        self.B_vals += self.B_vph.vals
        del self.B_vph

        self.add_gov_equation('B_thlorentz', 'bth')
        self.B_thlorentz.add_term('bth', ones, k_vals=range(1,Nk-1))
        self.B_rows += self.B_thlorentz.rows
        self.B_cols += self.B_thlorentz.cols
        self.B_vals += self.B_thlorentz.vals
        del self.B_thlorentz

        self.add_gov_equation('B_phlorentz', 'bph')
        self.B_phlorentz.add_term('bph', ones, k_vals=range(1,Nk-1))
        self.B_rows += self.B_phlorentz.rows
        self.B_cols += self.B_phlorentz.cols
        self.B_vals += self.B_phlorentz.vals
        del self.B_phlorentz

        self.add_gov_equation('B_ur', 'ur')
        self.B_ur.add_term('ur', ones)
        self.B_rows += self.B_ur.rows
        self.B_cols += self.B_ur.cols
        self.B_vals += self.B_ur.vals
        del self.B_ur
        self.B = FVF_model_base.coo_matrix((self.B_vals, (self.B_rows, self.B_cols)),
                                   shape=(self.SizeM, self.SizeM))
        del self.B_vals, self.B_rows, self.B_cols
        return self.B