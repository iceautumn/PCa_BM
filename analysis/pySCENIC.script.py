import argparse
import os
import loompy as lp
import pandas as pd

def scenic_workflow(loom, out, tfs_path, db_feather, motif_path, loom_output, AUCell_matrix_out):

    os.chdir(out)

    num_workers = 20

    outpath_adj = os.path.join(out, "adj.csv")

    os.system(
        f"pyscenic grn {loom} {tfs_path} -o {outpath_adj} --num_workers {num_workers}")

    os.system(
        f"pyscenic ctx {outpath_adj} {db_feather} --annotations_fname {motif_path} --expression_mtx_fname {loom} --output reg.csv --mask_dropouts --num_workers {num_workers} > pyscenic_ctx_stdout.txt")

    os.system(
        f"pyscenic aucell {loom} reg.csv --output {loom_output} --num_workers {num_workers} > pyscenic_aucell_stdout.txt")

    lf = lp.connect(loom_output, mode="r+", validate=False)
    auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    lf.close()

    auc_mtx.to_csv(AUCell_matrix_out, index=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Run pySCENIC workflow with input parameters.')
    parser.add_argument('--loom', type=str, help='The input loom file.')
    parser.add_argument('--out', type=str, help='The output directory.')
    parser.add_argument('--tfs_path', type=str, help='The tfs_path in the code.')
    parser.add_argument('--db_feather', type=str, help='The db_names in the code.')
    parser.add_argument('--motif_path', type=str, help='The motif_path for --annotations_fname.')
    parser.add_argument('--loom_output', type=str, default='processed_output.loom',
                        help='The loom_path_output in the code.')
    parser.add_argument('--AUCell_matrix_out', type=str, default='AUC_results.csv',
                        help='The output file for AUC results.')

    args = parser.parse_args()

    scenic_workflow(args.loom, args.out, args.tfs_path, args.db_feather, args.motif_path, args.loom_output, args.AUCell_matrix_out)

