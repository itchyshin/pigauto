# Data notes

Both datasets needed by the benchmark scripts are **bundled with the pigauto R
package** and do not require any manual upload or rsync.

Once `pigauto` is installed on Compute Canada (see the top-level README for
install instructions), the data objects are available via `utils::data()`:

```r
e <- new.env(parent = emptyenv())
utils::data("avonet_full", package = "pigauto", envir = e)
utils::data("tree_full",   package = "pigauto", envir = e)
```

Both R scripts in `Rscripts/` already contain these calls — no action needed.

## Dataset details

| Object       | Dimensions               | Format                    | Source                     |
|--------------|--------------------------|---------------------------|----------------------------|
| `avonet_full`| 9,993 species x 8 cols   | data.frame (includes `Species_Key`) | AVONET3 + BirdTree         |
| `tree_full`  | 9,993 tips               | `phylo` (ape)             | BirdTree maximum-clade-credibility tree |

The scripts align the dataset to the tree via `rownames(df) <- df$Species_Key`.

## If you need to add other data

If a future benchmark requires data not bundled in the package, copy the `.rds`
or `.csv` files into this `data/` directory and upload the whole
`submit_v090_cloud/` bundle to CC via the `rsync` command in the README.
