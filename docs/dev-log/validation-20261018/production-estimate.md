# Measured pilot: provisional production estimate

Planned jobs: 32240. Nominal workers: 48. Measured serial-equivalent work: 265.9 CPU-hours.
Provisional elapsed range: 5.5 to 14.8 hours at 48 single-thread workers, before queue time.
Timing-only inputs retain separate raw source identities. Successful measured routes supersede failed adapter startup timings; all timing attempts remain in pilot_timings.csv.
The upper planning estimate applies a 2x timing margin and 75% utilization. Runtime remains uncertain across scenarios and concurrent load. Scenarios without exact size matches use conservative cubic scaling. Failed/unavailable pilot routes remain flagged; their timings do not establish successful production feasibility.
Scenarios without a successful measured pilot: 2. Scenarios retaining a failed/unavailable timing attempt: 7. Size-extrapolated scenarios: 0. Mechanism-pooled timings: 25.
No production approval receipt has been created. Recheck host load, aggregate user CPU use and measured memory before scheduling.
