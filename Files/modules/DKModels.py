import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import *

class HomogeneousTransformationMatrix:
    
    def __init__(self, idx, wrt_idx) -> None:
        self.set_idx(idx)
        self.set_wrt_idx(wrt_idx)
        
        
    def set_idx(self, idx):
        self.idx = idx
    def set_wrt_idx(self, wrt_idx):
        self.wrt_idx = wrt_idx    

class DHframe:
    def __init__(self, idx, wrt_idx) -> None:
        
        self.set_idx(idx)
        self.set_wrt_idx(wrt_idx)
        
        self.d, self.a, self.teta, self.alpha = symbols(f'd{idx}, a{idx}, theta{idx},alpha{idx}')

        self.create_DH_matrix()


    def create_DH_matrix(self):
        d, a, teta, alpha = self.d, self.a, self.teta, self.alpha
        DH_matrix = Matrix([
            [cos(teta), -cos(alpha)*sin(teta), sin(alpha)*sin(teta), a*cos(teta)],
            [sin(teta), cos(alpha)*cos(teta), -sin(alpha)*cos(teta), a*sin(teta)],
            [0, sin(alpha), cos(alpha), d],
            [0, 0, 0, 1]
        ])
        self.DH_matrix = DH_matrix
            
    
    def change_d(self, new_d):
        self.DH_matrix = self.DH_matrix.subs(self.d, new_d)
        self.d = new_d
    
    def change_a(self, new_a):
        self.DH_matrix = self.DH_matrix.subs(self.a, new_a)
        self.a = new_a
        
    def change_teta(self, new_teta):
        self.DH_matrix = self.DH_matrix.subs(self.teta, new_teta)
        self.teta = new_teta
        
    def change_alpha(self, new_alpha):
        self.DH_matrix = self.DH_matrix.subs(self.alpha, new_alpha)
        self.alpha = new_alpha
    
    def fill_DH_matrix(self, d, a, teta, alpha):
        self.change_d(d)
        self.change_a(a)
        self.change_teta(teta)
        self.change_alpha(alpha)


    def set_idx(self, idx):
        self.idx = idx
    def set_wrt_idx(self, wrt_idx):
        self.wrt_idx = wrt_idx
    
    def get_DH_matrix(self):
        return self.DH_matrix
    
    def get_position(self):
        return self.DH_matrix[0:3, 3]
    
    def get_orientation(self):
        return self.DH_matrix[0:3, 0:3]
        
        
class DireckKinematic:
    
    def __init__(self, DH_table) -> None:
        
        self.DH_table = DH_table
        self.DH_frames = []
        self.create_frames()
        
    
    def add_frame(self, frame, idx):
        self.DH_frames.insert(idx, frame)
    
    def remove_frame(self, idx):
        self.DH_frames.pop(idx)

    def create_frames(self):
        for i, params in enumerate(self.DH_table):
            d, a, teta, alpha= params
            frame = DHframe(idx=i+1, wrt_idx=i)
            frame.fill_DH_matrix(d, a, teta, alpha)
            self.add_frame(frame, i)
            
        
    
    def compute_DK(self, base_idx, end_idx):
        frames = self.DH_frames
        mat = frames[base_idx].get_DH_matrix()
        for i in range(base_idx+1, end_idx):
            mat = mat * frames[i].get_DH_matrix()
        
        mat = trigsimp(mat)
        rotation = mat[0:3, 0:3]
        position = mat[0:3, 3]
        return mat, rotation, position
        
    
    def apply_configuration(self, config, base_idx, end_idx):
        mat, rotation, position = self.compute_DK(base_idx, end_idx)
        mat = mat.subs(config)
        rotation = rotation.subs(config)
        position = position.subs(config)
        return mat, rotation, position
        
    def show_DH_table(self):
        dh_tb = pd.DataFrame(self.DH_table, columns=['d', 'a', 'theta', 'alpha'])
        print(dh_tb)