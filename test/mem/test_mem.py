import os
import shutil

# def test_trjcat(tmp_path):
#     pathpath= os.path.dirname(__file__) + "/../simdata_dw/load/210"
#
# def test_trjcat(tmp_path):
#     pathpath= os.path.dirname(__file__) + "/../simdata_dw/load/210"

def test_mem_com(tmp_path):
    from dztools.mem.cvs import mem_helicity
    mlt_f = os.path.dirname(__file__) + "/data/2mlt.pdb"
    mlt_u = os.path.dirname(__file__) + "/data/2mlt_u.pdb"

    hel_f = mem_helicity(top=mlt_f, xtc=mlt_f)
    hel_u = mem_helicity(top=mlt_u, xtc=mlt_u)

    assert abs(hel_f[0]-0.92) < 0.0001
    assert abs(hel_u[0]-0.00) < 0.0001

def test_mem_helicity(tmp_path):
    from dztools.mem.cvs import mem_helicity
    mlt_f = os.path.dirname(__file__) + "/data/2mlt.pdb"
    mlt_u = os.path.dirname(__file__) + "/data/2mlt_u.pdb"

    hel_f = mem_helicity(top=mlt_f, xtc=mlt_f)
    hel_u = mem_helicity(top=mlt_u, xtc=mlt_u)

    assert abs(hel_f[0]-0.92) < 0.0001
    assert abs(hel_u[0]-0.00) < 0.0001

def test_mem_chain(tmp_path):
    from dztools.mem.cvs import mem_chain
    pull_s = os.path.dirname(__file__) + "/data/pull_s.gro"
    pull_e = os.path.dirname(__file__) + "/data/pull_e.gro"

    hoxy = "OW OA OB OC OD"
    lip = "DMPC"

    epsilon_s = mem_chain(top=pull_s, xtc=pull_s, hoxy=hoxy, lip=lip)
    epsilon_e = mem_chain(top=pull_e, xtc=pull_e, hoxy=hoxy, lip=lip)

    assert abs(epsilon_s-0.1901430306128322) < 0.0001
    assert abs(epsilon_e-0.7483991044081285) < 0.0001

