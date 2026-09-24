# Four-tone BER qualification results

Both frozen channels pass the post-decoding BER target of 1e-3 with joint coverage of at least 95%, using time-uniform packet-fraction confidence sequences and alpha=0.025 per channel. All tentative payload errors, including CRC-rejected packets, are counted. Screening packets are excluded.

| Links | Carrier | Packets | Payload errors / bits | BER | Confidence interval | CRC rejections | Wrong accepted | Max latency |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| 0.5 m and 1 m | 24 Hz | 16,896 | 2,303 / 3,784,704 | 0.000608502 | [0.000437218, 0.000997235] | 237 | 0 | 224.098 s |
| 2 m | 18 Hz | 9,344 | 666 / 2,093,056 | 0.000318195 | [0.000167922, 0.000989752] | 73 | 0 | 224.100 s |

Nominal payload rate is 1 bit/s. CRC-accepted payload rates are 0.98597 bit/s (24 Hz) and 0.99219 bit/s (18 Hz). Maximum measured processing blocks are 0.1801 s and 0.1898 s respectively. Timing measurements include parallel host load and do not certify deployment worst-case performance.

The current rotor sizes, source amplitudes and receiver resources were held fixed. Four tones reduce peak inertial modulation torque by about 16% relative to the qualified seven-tone protocol at these same rotor sizes. This study does not establish minimum transmitter size.

Verification: source and snapshot hashes match the frozen manifest; packet indices form complete independent-seed prefixes; final statistics were recomputed from all 26,240 packet records; every earlier 128-packet stopping prefix remains inconclusive. All per-packet detector quality checks and two fresh wider-history/doubled-lookahead convergence checks per channel passed.

Outstanding detector assessment: wider-detector replays of qualification error packets have not been performed. Their seeds are saved in verification.json. Statistical qualification applies to the implemented finite-history streaming receiver; two diagnostic packets per channel do not establish convergence for every error event. Acquisition is budgeted, with initial timing, phase and state known in the simulation. Zero wrong accepted packets is an observation, not a zero undetected-error guarantee.

The manuscript remains unchanged pending review.
