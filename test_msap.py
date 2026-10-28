from msap_sim import candidates,rc,cut_sites

def test_fragment_orientation():
    # E left, H right; E extension AAC; H extension = RC(TGA) = TCA
    seq='GAATTC'+'AAC'+'A'*150+'TGA'+'CCGG'
    f=list(candidates('test',seq,100,500,3))
    assert len(f)==1 and f[0]['E_extension']=='AAC' and f[0]['HM_extension']=='TCA',f
    # H left, E right; H extension GGT; E extension RC(CCA) = TGG
    seq='CCGG'+'GGT'+'A'*150+'CCA'+'GAATTC'
    f=list(candidates('test',seq,100,500,3))
    assert len(f)==1 and f[0]['HM_extension']=='GGT' and f[0]['E_extension']=='TGG',f

def test_internal_cut_breaks_fragment():
    seq='GAATTC'+'A'*150+'CCGG'+'A'*140+'CCGG'
    f=list(candidates('test',seq,100,500,3))
    assert len(f)==1,f

def test_rc(): assert rc('ATCG')=='CGAT'

if __name__=='__main__':
    test_fragment_orientation();test_internal_cut_breaks_fragment();test_rc(); print('3 tests passed')
