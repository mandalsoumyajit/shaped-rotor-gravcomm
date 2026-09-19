import numpy as np
from gravcomm.packet_confidence import packet_ber_interval


def test_zero_and_all_error_limits():
    zero=packet_ber_interval(np.zeros(1000))
    one=packet_ber_interval(np.full(1000,64))
    assert zero[0]==0 and 0<zero[1]<.01
    assert one[1]==1
    np.testing.assert_allclose(one,[1-zero[1],1],rtol=1e-10)


def test_burst_correlation_and_order_invariance():
    burst=np.r_[np.zeros(990),np.full(10,8)]
    interval=packet_ber_interval(burst)
    assert interval[0]<burst.mean()/64<interval[1]
    np.testing.assert_allclose(interval,packet_ber_interval(burst[::-1]))
    # Holding the number of bit errors fixed while spreading them among more
    # independent packets carries more information about the mean.
    spread=np.r_[np.zeros(920),np.ones(80)]
    better=packet_ber_interval(spread)
    assert better[0]>interval[0]
