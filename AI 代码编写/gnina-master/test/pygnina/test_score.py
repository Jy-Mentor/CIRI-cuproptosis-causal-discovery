import pytest
import pygnina

from data import *

def test_score():
  G = pygnina.GNINA()
  G.set_receptor(REC,'pdb')
  R = G.score(LIGSDF,'sdf')
  assert R.energy() == pytest.approx(-8.2394, rel=1e-3)
  assert R.cnnscore() == pytest.approx(0.9813, rel=1e-3)
  assert R.cnnaffinity() == pytest.approx(7.247, rel=1e-3)






