# scPCIP container image

One Apptainer / Singularity image with everything `scpcip` needs, including the optional R figures
(organ maps drawn with ggPlantmap, UMAPs and marker heatmaps drawn with SCP) and the gene -> orthogroup map,
so it runs offline.

| Component | Version |
|---|---|
| Python | 3.13 |
| R | 4.5.3 |
| Seurat / SeuratObject | 5.5.0 / 5.4.0 |
| SCP | GitHub zhanghao-njmu/SCP @ b9b0eb7 (0.5.6) |
| ggPlantmap | GitHub leonardojo/ggPlantmap @ e050817 |
| gene -> orthogroup map | Hugging Face revision 7eaec1bb7081 |

## Use

```bash
apptainer run scpcip.sif predict --input query.h5ad --species "Zea mays" --out results/
apptainer run scpcip.sif plot --input query.h5ad --species "Zea mays" --out results/
```

`apptainer run scpcip.sif <args>` is `scpcip <args>`. Your home and current directory are mounted by default;
mount other data directories with `--bind /path/to/data`.

## Build

```bash
mkdir build && cd build
cp ../container/scpcip.def .
rsync -a --exclude build --exclude examples --exclude container ../ scpcip/
scpcip fetch                                         # prints the path of the cached map
cp "$(scpcip fetch | tail -1)" gene_to_feature.tsv
mkdir rpkgs
curl -fsSL -o rpkgs/SCP_b9b0eb7a7bf2c2c4b2262e73e09d7ebd515c7da0.tar.gz \
  https://codeload.github.com/zhanghao-njmu/SCP/tar.gz/b9b0eb7a7bf2c2c4b2262e73e09d7ebd515c7da0
curl -fsSL -o rpkgs/ggPlantmap_e0508172709dc192dfc3e99da631984dbd9a97d2.tar.gz \
  https://codeload.github.com/leonardojo/ggPlantmap/tar.gz/e0508172709dc192dfc3e99da631984dbd9a97d2
mkdir fonts
curl -fsSL https://github.com/dejavu-fonts/dejavu-fonts/releases/download/version_2_37/dejavu-fonts-ttf-2.37.tar.bz2 \
  | tar xj --strip-components=2 -C fonts dejavu-fonts-ttf-2.37/ttf/DejaVuSans.ttf dejavu-fonts-ttf-2.37/ttf/DejaVuSans-Bold.ttf
apptainer build scpcip.sif scpcip.def                # add --fakeroot when not root
```

Without root, if `--fakeroot` stops with "newuidmap was not found", build inside a user namespace instead:
`unshare -r apptainer build scpcip.sif scpcip.def`.
