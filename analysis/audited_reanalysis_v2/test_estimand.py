"""Small mathematical QA fixtures only; never used as study observations."""
import numpy as np
from analyze_verified import aipw

def test_exact_ipw_reduction():
    a=np.array([0.,0.,1.,1.]);y=np.array([0.,1.,0.,1.]);e=np.array([.25,.25,.25,.25]);m=np.zeros(4)
    h1,h0=aipw(a,y,e,m,m)
    assert np.allclose(h1,a*y/e)
    assert np.allclose(h0,(1-a)*y/(1-e))
    assert np.isclose(h1.mean()-h0.mean(),2/3)

def test_zero_residual_outcome_model():
    a=np.array([0.,1.]);y=np.array([.2,.7]);e=np.array([.2,.8])
    h1,h0=aipw(a,y,e,np.array([.7,.7]),np.array([.2,.2]))
    assert np.isclose((h1-h0).mean(),.5)

if __name__=='__main__':
    test_exact_ipw_reduction();test_zero_residual_outcome_model();print('AIPW formula QA passed (2 fixtures)')
