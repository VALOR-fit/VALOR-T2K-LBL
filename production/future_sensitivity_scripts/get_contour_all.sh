# Loop over all command line inputs (all assumed to be file names)
for f in "$@"; do
    fit_path=$(dirname "${f}")
    tag=$(basename "${f}")
    tag="${tag#*fit}"
    tag="${tag%-*}"
    tag="get_contour${tag}"
    ${T2KOSCANA}/production/submit_get_contour --file ${f} --resource interactive --job-dir ${fit_path} --tag ${tag}
done
