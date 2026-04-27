## Publication

ArXiV: 

## Tips & Tricks

1. If running this on your computing cluster, consider making a `model.stan` and a separate `model-nogq.stan`, where the latter has no generated quantities block. The GQs for these models are very expensive to store. With the two files separate, you can run `model-gq.stan` first and run the generated quantities block of `model.stan` separately-- this is why [script.R] is written to delete GQ draws immediately after summarizing them. 