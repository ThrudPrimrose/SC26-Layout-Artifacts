# Array co-access analysis: `velocity_tendencies`

34 unique arrays across 17 loop nests, forming 72 array-nest pairs.

## 1. Per-nest working sets

Each row lists the arrays that a single loop nest accesses together. Arrays within the same row share a working set: they are loaded into cache during the same inner-loop execution and therefore benefit from co-located storage.

| Nest | Shape | Collapsed | Behavior | Size | Arrays |
|-----:|-------|-----------|----------|-----:|--------|
| 1 | `b.h` S | `h:S  b:S` |  | 6 | `p_diag%vn_ie`, `p_diag%vt`, `p_metrics%wgtfacq_e`, `p_prog%vn`, `z_kin_hor_e`, `z_vt_ie` |
| 2 | `b.h` S | `h:S  b:S` | accumulate | 7 | `p_diag%vn_ie`, `p_diag%vn_ie_ubc`, `p_diag%vt`, `p_metrics%wgtfacq_e`, `p_prog%vn`, `z_kin_hor_e`, `z_vt_ie` |
| 3 | `v.h` **U** | `h:U  v:S  b:S` |  | 3 | `p_diag%vt`, `p_int%rbf_vec_coeff_e`, `p_prog%vn` |
| 4 | `v.h` S | `h:S  v:S  b:S` |  | 5 | `p_diag%vn_ie`, `p_diag%vt`, `p_metrics%wgtfac_e`, `p_prog%vn`, `z_kin_hor_e` |
| 5 | `v.h` S | `h:S  v:S  b:S` |  | 3 | `p_diag%vt`, `p_metrics%wgtfac_e`, `z_vt_ie` |
| 6 | `v.h` S | `h:S  v:S  b:S` |  | 5 | `p_diag%vt`, `p_metrics%ddxn_z_full`, `p_metrics%ddxt_z_full`, `p_prog%vn`, `z_w_concorr_me` |
| 7 | `v.h` **U** | `h:U  v:S  b:S` |  | 8 | `p_diag%vn_ie`, `p_patch%edges%inv_dual_edge_length`, `p_patch%edges%inv_primal_edge_length`, `p_patch%edges%tangent_orientation`, `p_prog%w`, `z_v_grad_w`, `z_vt_ie`, `z_w_v` |
| 8 | `b.h` S | `h:S` |  | 1 | `z_w_con_c` |
| 9 | `v.h` **U** | `h:U  v:S  b:S` |  | 3 | `p_int%e_bln_c_s`, `z_ekinh`, `z_kin_hor_e` |
| 10 | `v.h` **U** | `h:U  v:S  b:S` |  | 3 | `p_int%e_bln_c_s`, `z_w_concorr_mc`, `z_w_concorr_me` |
| 11 | `v.h` S | `h:S  v:S  b:S` |  | 3 | `p_diag%w_concorr_c`, `p_metrics%wgtfac_c`, `z_w_concorr_mc` |
| 12 | `v.h` S | `h:S  v:S  b:S` |  | 2 | `p_prog%w`, `z_w_con_c` |
| 13 | `v.h` S | `h:S  v:S  b:S` | accumulate | 2 | `p_diag%w_concorr_c`, `z_w_con_c` |
| 14 | `v.h` S | `h:S  v:S  b:S` |  | 2 | `z_w_con_c`, `z_w_con_c_full` |
| 15 | `v.h` S | `h:S  v:S  b:S` |  | 5 | `p_diag%ddt_w_adv_pc`, `p_metrics%coeff1_dwdz`, `p_metrics%coeff2_dwdz`, `p_prog%w`, `z_w_con_c` |
| 16 | `v.h` **U** | `h:U  v:S  b:S` | accumulate | 3 | `p_diag%ddt_w_adv_pc`, `p_int%e_bln_c_s`, `z_v_grad_w` |
| 17 | `v.h` **U** | `h:U  v:S  b:S` |  | 11 | `p_diag%ddt_vn_apc_pc`, `p_diag%vn_ie`, `p_diag%vt`, `p_int%c_lin_e`, `p_metrics%coeff_gradekin`, `p_metrics%ddqz_z_full_e`, `p_patch%edges%f_e`, `z_ekinh`, `z_kin_hor_e`, `z_w_con_c_full`, `zeta` |

## 2. Per-array nest membership

Each row shows which nests access a given array. Arrays appearing in many nests have the most constrained layout requirements.

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%vt` | 7 | [1, 2, 3, 4, 5, 6, 17] |
| `p_diag%vn_ie` | 5 | [1, 2, 4, 7, 17] |
| `p_prog%vn` | 5 | [1, 2, 3, 4, 6] |
| `z_kin_hor_e` | 5 | [1, 2, 4, 9, 17] |
| `z_w_con_c` | 5 | [8, 12, 13, 14, 15] |
| `z_vt_ie` | 4 | [1, 2, 5, 7] |
| `p_int%e_bln_c_s` | 3 | [9, 10, 16] |
| `p_prog%w` | 3 | [7, 12, 15] |
| `p_diag%ddt_w_adv_pc` | 2 | [15, 16] |
| `p_diag%w_concorr_c` | 2 | [11, 13] |
| `p_metrics%wgtfac_e` | 2 | [4, 5] |
| `p_metrics%wgtfacq_e` | 2 | [1, 2] |
| `z_ekinh` | 2 | [9, 17] |
| `z_v_grad_w` | 2 | [7, 16] |
| `z_w_con_c_full` | 2 | [14, 17] |
| `z_w_concorr_mc` | 2 | [10, 11] |
| `z_w_concorr_me` | 2 | [6, 10] |
| `p_diag%ddt_vn_apc_pc` | 1 | [17] |
| `p_diag%vn_ie_ubc` | 1 | [2] |
| `p_int%c_lin_e` | 1 | [17] |
| `p_int%rbf_vec_coeff_e` | 1 | [3] |
| `p_metrics%coeff1_dwdz` | 1 | [15] |
| `p_metrics%coeff2_dwdz` | 1 | [15] |
| `p_metrics%coeff_gradekin` | 1 | [17] |
| `p_metrics%ddqz_z_full_e` | 1 | [17] |
| `p_metrics%ddxn_z_full` | 1 | [6] |
| `p_metrics%ddxt_z_full` | 1 | [6] |
| `p_metrics%wgtfac_c` | 1 | [11] |
| `p_patch%edges%f_e` | 1 | [17] |
| `p_patch%edges%inv_dual_edge_length` | 1 | [7] |
| `p_patch%edges%inv_primal_edge_length` | 1 | [7] |
| `p_patch%edges%tangent_orientation` | 1 | [7] |
| `z_w_v` | 1 | [7] |
| `zeta` | 1 | [17] |

## 3. Pairwise co-occurrence

Of the 561 possible array pairs, 139 share at least one nest.

| Array A | Array B | Shared nests | Which nests |
|---------|---------|:------------:|-------------|
| `p_diag%vt` | `p_prog%vn` | 5 | [1, 2, 3, 4, 6] |
| `p_diag%vn_ie` | `p_diag%vt` | 4 | [1, 2, 4, 17] |
| `p_diag%vn_ie` | `z_kin_hor_e` | 4 | [1, 2, 4, 17] |
| `p_diag%vt` | `z_kin_hor_e` | 4 | [1, 2, 4, 17] |
| `p_diag%vn_ie` | `p_prog%vn` | 3 | [1, 2, 4] |
| `p_diag%vn_ie` | `z_vt_ie` | 3 | [1, 2, 7] |
| `p_diag%vt` | `z_vt_ie` | 3 | [1, 2, 5] |
| `p_prog%vn` | `z_kin_hor_e` | 3 | [1, 2, 4] |
| `p_diag%vn_ie` | `p_metrics%wgtfacq_e` | 2 | [1, 2] |
| `p_diag%vt` | `p_metrics%wgtfac_e` | 2 | [4, 5] |
| `p_diag%vt` | `p_metrics%wgtfacq_e` | 2 | [1, 2] |
| `p_metrics%wgtfacq_e` | `p_prog%vn` | 2 | [1, 2] |
| `p_metrics%wgtfacq_e` | `z_kin_hor_e` | 2 | [1, 2] |
| `p_metrics%wgtfacq_e` | `z_vt_ie` | 2 | [1, 2] |
| `p_prog%vn` | `z_vt_ie` | 2 | [1, 2] |
| `p_prog%w` | `z_w_con_c` | 2 | [12, 15] |
| `z_ekinh` | `z_kin_hor_e` | 2 | [9, 17] |
| `z_kin_hor_e` | `z_vt_ie` | 2 | [1, 2] |
| `p_diag%ddt_vn_apc_pc` | `p_diag%vn_ie` | 1 | [17] |
| `p_diag%ddt_vn_apc_pc` | `p_diag%vt` | 1 | [17] |
| `p_diag%ddt_vn_apc_pc` | `p_int%c_lin_e` | 1 | [17] |
| `p_diag%ddt_vn_apc_pc` | `p_metrics%coeff_gradekin` | 1 | [17] |
| `p_diag%ddt_vn_apc_pc` | `p_metrics%ddqz_z_full_e` | 1 | [17] |
| `p_diag%ddt_vn_apc_pc` | `p_patch%edges%f_e` | 1 | [17] |
| `p_diag%ddt_vn_apc_pc` | `z_ekinh` | 1 | [17] |
| `p_diag%ddt_vn_apc_pc` | `z_kin_hor_e` | 1 | [17] |
| `p_diag%ddt_vn_apc_pc` | `z_w_con_c_full` | 1 | [17] |
| `p_diag%ddt_vn_apc_pc` | `zeta` | 1 | [17] |
| `p_diag%ddt_w_adv_pc` | `p_int%e_bln_c_s` | 1 | [16] |
| `p_diag%ddt_w_adv_pc` | `p_metrics%coeff1_dwdz` | 1 | [15] |
| `p_diag%ddt_w_adv_pc` | `p_metrics%coeff2_dwdz` | 1 | [15] |
| `p_diag%ddt_w_adv_pc` | `p_prog%w` | 1 | [15] |
| `p_diag%ddt_w_adv_pc` | `z_v_grad_w` | 1 | [16] |
| `p_diag%ddt_w_adv_pc` | `z_w_con_c` | 1 | [15] |
| `p_diag%vn_ie` | `p_diag%vn_ie_ubc` | 1 | [2] |
| `p_diag%vn_ie` | `p_int%c_lin_e` | 1 | [17] |
| `p_diag%vn_ie` | `p_metrics%coeff_gradekin` | 1 | [17] |
| `p_diag%vn_ie` | `p_metrics%ddqz_z_full_e` | 1 | [17] |
| `p_diag%vn_ie` | `p_metrics%wgtfac_e` | 1 | [4] |
| `p_diag%vn_ie` | `p_patch%edges%f_e` | 1 | [17] |
| `p_diag%vn_ie` | `p_patch%edges%inv_dual_edge_length` | 1 | [7] |
| `p_diag%vn_ie` | `p_patch%edges%inv_primal_edge_length` | 1 | [7] |
| `p_diag%vn_ie` | `p_patch%edges%tangent_orientation` | 1 | [7] |
| `p_diag%vn_ie` | `p_prog%w` | 1 | [7] |
| `p_diag%vn_ie` | `z_ekinh` | 1 | [17] |
| `p_diag%vn_ie` | `z_v_grad_w` | 1 | [7] |
| `p_diag%vn_ie` | `z_w_con_c_full` | 1 | [17] |
| `p_diag%vn_ie` | `z_w_v` | 1 | [7] |
| `p_diag%vn_ie` | `zeta` | 1 | [17] |
| `p_diag%vn_ie_ubc` | `p_diag%vt` | 1 | [2] |

(89 additional pairs omitted.)

## 4. Co-occurrence frequency distribution

How many array pairs share exactly k nests.

| Shared nests (k) | #Pairs | Fraction |
|-----------------:|-------:|---------:|
| 0 | 422 | 75.2% |
| 1 | 121 | 21.6% |
| 2 | 10 | 1.8% |
| 3 | 4 | 0.7% |
| 4 | 3 | 0.5% |
| 5 | 1 | 0.2% |

## 5. Equivalence classes (identical nest membership)

Arrays that appear in exactly the same set of nests. They are indistinguishable from a co-access perspective.

### Multi-member classes (4 classes, 14 arrays)

**Nests {17}** — 6 arrays:

- `p_diag%ddt_vn_apc_pc`
- `p_int%c_lin_e`
- `p_metrics%coeff_gradekin`
- `p_metrics%ddqz_z_full_e`
- `p_patch%edges%f_e`
- `zeta`

**Nests {7}** — 4 arrays:

- `p_patch%edges%inv_dual_edge_length`
- `p_patch%edges%inv_primal_edge_length`
- `p_patch%edges%tangent_orientation`
- `z_w_v`

**Nests {6}** — 2 arrays:

- `p_metrics%ddxn_z_full`
- `p_metrics%ddxt_z_full`

**Nests {15}** — 2 arrays:

- `p_metrics%coeff1_dwdz`
- `p_metrics%coeff2_dwdz`

### Singletons (20 arrays)

| Array | Nests |
|-------|-------|
| `p_metrics%wgtfacq_e` | [1, 2] |
| `p_diag%vt` | [1, 2, 3, 4, 5, 6, 17] |
| `p_prog%vn` | [1, 2, 3, 4, 6] |
| `p_diag%vn_ie` | [1, 2, 4, 7, 17] |
| `z_kin_hor_e` | [1, 2, 4, 9, 17] |
| `z_vt_ie` | [1, 2, 5, 7] |
| `p_diag%vn_ie_ubc` | [2] |
| `p_int%rbf_vec_coeff_e` | [3] |
| `p_metrics%wgtfac_e` | [4, 5] |
| `z_w_concorr_me` | [6, 10] |
| `p_prog%w` | [7, 12, 15] |
| `z_v_grad_w` | [7, 16] |
| `z_w_con_c` | [8, 12, 13, 14, 15] |
| `p_int%e_bln_c_s` | [9, 10, 16] |
| `z_ekinh` | [9, 17] |
| `z_w_concorr_mc` | [10, 11] |
| `p_metrics%wgtfac_c` | [11] |
| `p_diag%w_concorr_c` | [11, 13] |
| `z_w_con_c_full` | [14, 17] |
| `p_diag%ddt_w_adv_pc` | [15, 16] |

## 6. Nest intersection sizes

How many arrays two nests share.

| | N1 | N2 | N3 | N4 | N5 | N6 | N7 | N8 | N9 | N10 | N11 | N12 | N13 | N14 | N15 | N16 | N17 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **N1** | 6 | 6 | 2 | 4 | 2 | 2 | 2 | · | 1 | · | · | · | · | · | · | · | 3 |
| **N2** | 6 | 7 | 2 | 4 | 2 | 2 | 2 | · | 1 | · | · | · | · | · | · | · | 3 |
| **N3** | 2 | 2 | 3 | 2 | 1 | 2 | · | · | · | · | · | · | · | · | · | · | 1 |
| **N4** | 4 | 4 | 2 | 5 | 2 | 2 | 1 | · | 1 | · | · | · | · | · | · | · | 3 |
| **N5** | 2 | 2 | 1 | 2 | 3 | 1 | 1 | · | · | · | · | · | · | · | · | · | 1 |
| **N6** | 2 | 2 | 2 | 2 | 1 | 5 | · | · | · | 1 | · | · | · | · | · | · | 1 |
| **N7** | 2 | 2 | · | 1 | 1 | · | 8 | · | · | · | · | 1 | · | · | 1 | 1 | 1 |
| **N8** | · | · | · | · | · | · | · | 1 | · | · | · | 1 | 1 | 1 | 1 | · | · |
| **N9** | 1 | 1 | · | 1 | · | · | · | · | 3 | 1 | · | · | · | · | · | 1 | 2 |
| **N10** | · | · | · | · | · | 1 | · | · | 1 | 3 | 1 | · | · | · | · | 1 | · |
| **N11** | · | · | · | · | · | · | · | · | · | 1 | 3 | · | 1 | · | · | · | · |
| **N12** | · | · | · | · | · | · | 1 | 1 | · | · | · | 2 | 1 | 1 | 2 | · | · |
| **N13** | · | · | · | · | · | · | · | 1 | · | · | 1 | 1 | 2 | 1 | 1 | · | · |
| **N14** | · | · | · | · | · | · | · | 1 | · | · | · | 1 | 1 | 2 | 1 | · | 1 |
| **N15** | · | · | · | · | · | · | 1 | 1 | · | · | · | 2 | 1 | 1 | 5 | 1 | · |
| **N16** | · | · | · | · | · | · | 1 | · | 1 | 1 | · | · | · | · | 1 | 3 | · |
| **N17** | 3 | 3 | 1 | 3 | 1 | 1 | 1 | · | 2 | · | · | · | · | 1 | · | · | 11 |

## 7. Jaccard similarity (top pairs)

Jaccard index J(A,B) = |nests(A) ∩ nests(B)| / |nests(A) ∪ nests(B)| normalises co-occurrence by the combined footprint of both arrays. J = 1 means identical nest membership; J = 0 means disjoint.

Of 561 pairs: 52 have J >= 0.5, 422 have J = 0 (disjoint).

| Array A | Array B | Jaccard | Shared | |A| | |B| |
|---------|---------|:-------:|:------:|----:|----:|
| `p_diag%ddt_vn_apc_pc` | `p_int%c_lin_e` | 1.00 | 1 | 1 | 1 |
| `p_diag%ddt_vn_apc_pc` | `p_metrics%coeff_gradekin` | 1.00 | 1 | 1 | 1 |
| `p_diag%ddt_vn_apc_pc` | `p_metrics%ddqz_z_full_e` | 1.00 | 1 | 1 | 1 |
| `p_diag%ddt_vn_apc_pc` | `p_patch%edges%f_e` | 1.00 | 1 | 1 | 1 |
| `p_diag%ddt_vn_apc_pc` | `zeta` | 1.00 | 1 | 1 | 1 |
| `p_int%c_lin_e` | `p_metrics%coeff_gradekin` | 1.00 | 1 | 1 | 1 |
| `p_int%c_lin_e` | `p_metrics%ddqz_z_full_e` | 1.00 | 1 | 1 | 1 |
| `p_int%c_lin_e` | `p_patch%edges%f_e` | 1.00 | 1 | 1 | 1 |
| `p_int%c_lin_e` | `zeta` | 1.00 | 1 | 1 | 1 |
| `p_metrics%coeff1_dwdz` | `p_metrics%coeff2_dwdz` | 1.00 | 1 | 1 | 1 |
| `p_metrics%coeff_gradekin` | `p_metrics%ddqz_z_full_e` | 1.00 | 1 | 1 | 1 |
| `p_metrics%coeff_gradekin` | `p_patch%edges%f_e` | 1.00 | 1 | 1 | 1 |
| `p_metrics%coeff_gradekin` | `zeta` | 1.00 | 1 | 1 | 1 |
| `p_metrics%ddqz_z_full_e` | `p_patch%edges%f_e` | 1.00 | 1 | 1 | 1 |
| `p_metrics%ddqz_z_full_e` | `zeta` | 1.00 | 1 | 1 | 1 |
| `p_metrics%ddxn_z_full` | `p_metrics%ddxt_z_full` | 1.00 | 1 | 1 | 1 |
| `p_patch%edges%f_e` | `zeta` | 1.00 | 1 | 1 | 1 |
| `p_patch%edges%inv_dual_edge_length` | `p_patch%edges%inv_primal_edge_length` | 1.00 | 1 | 1 | 1 |
| `p_patch%edges%inv_dual_edge_length` | `p_patch%edges%tangent_orientation` | 1.00 | 1 | 1 | 1 |
| `p_patch%edges%inv_dual_edge_length` | `z_w_v` | 1.00 | 1 | 1 | 1 |
| `p_patch%edges%inv_primal_edge_length` | `p_patch%edges%tangent_orientation` | 1.00 | 1 | 1 | 1 |
| `p_patch%edges%inv_primal_edge_length` | `z_w_v` | 1.00 | 1 | 1 | 1 |
| `p_patch%edges%tangent_orientation` | `z_w_v` | 1.00 | 1 | 1 | 1 |
| `p_diag%vt` | `p_prog%vn` | 0.71 | 5 | 7 | 5 |
| `p_diag%vn_ie` | `z_kin_hor_e` | 0.67 | 4 | 5 | 5 |
| `p_diag%ddt_vn_apc_pc` | `z_ekinh` | 0.50 | 1 | 1 | 2 |
| `p_diag%ddt_vn_apc_pc` | `z_w_con_c_full` | 0.50 | 1 | 1 | 2 |
| `p_diag%ddt_w_adv_pc` | `p_metrics%coeff1_dwdz` | 0.50 | 1 | 2 | 1 |
| `p_diag%ddt_w_adv_pc` | `p_metrics%coeff2_dwdz` | 0.50 | 1 | 2 | 1 |
| `p_diag%vn_ie` | `p_diag%vt` | 0.50 | 4 | 5 | 7 |
| `p_diag%vn_ie` | `z_vt_ie` | 0.50 | 3 | 5 | 4 |
| `p_diag%vn_ie_ubc` | `p_metrics%wgtfacq_e` | 0.50 | 1 | 1 | 2 |
| `p_diag%vt` | `z_kin_hor_e` | 0.50 | 4 | 7 | 5 |
| `p_diag%w_concorr_c` | `p_metrics%wgtfac_c` | 0.50 | 1 | 2 | 1 |
| `p_int%c_lin_e` | `z_ekinh` | 0.50 | 1 | 1 | 2 |
| `p_int%c_lin_e` | `z_w_con_c_full` | 0.50 | 1 | 1 | 2 |
| `p_metrics%coeff_gradekin` | `z_ekinh` | 0.50 | 1 | 1 | 2 |
| `p_metrics%coeff_gradekin` | `z_w_con_c_full` | 0.50 | 1 | 1 | 2 |
| `p_metrics%ddqz_z_full_e` | `z_ekinh` | 0.50 | 1 | 1 | 2 |
| `p_metrics%ddqz_z_full_e` | `z_w_con_c_full` | 0.50 | 1 | 1 | 2 |
| `p_metrics%ddxn_z_full` | `z_w_concorr_me` | 0.50 | 1 | 1 | 2 |
| `p_metrics%ddxt_z_full` | `z_w_concorr_me` | 0.50 | 1 | 1 | 2 |
| `p_metrics%wgtfac_c` | `z_w_concorr_mc` | 0.50 | 1 | 1 | 2 |
| `p_metrics%wgtfacq_e` | `z_vt_ie` | 0.50 | 2 | 2 | 4 |
| `p_patch%edges%f_e` | `z_ekinh` | 0.50 | 1 | 1 | 2 |
| `p_patch%edges%f_e` | `z_w_con_c_full` | 0.50 | 1 | 1 | 2 |
| `p_patch%edges%inv_dual_edge_length` | `z_v_grad_w` | 0.50 | 1 | 1 | 2 |
| `p_patch%edges%inv_primal_edge_length` | `z_v_grad_w` | 0.50 | 1 | 1 | 2 |
| `p_patch%edges%tangent_orientation` | `z_v_grad_w` | 0.50 | 1 | 1 | 2 |
| `z_ekinh` | `zeta` | 0.50 | 1 | 2 | 1 |

### Jaccard distribution

| Jaccard range | #Pairs | Fraction |
|:-------------:|-------:|---------:|
| 0.00 | 422 | 75.2% |
| 0.01–0.24 | 53 | 9.4% |
| 0.25–0.49 | 34 | 6.1% |
| 0.50–0.74 | 29 | 5.2% |
| 0.75–0.99 | 0 | 0.0% |
| 1.00 | 23 | 4.1% |

## 8. Greedy Jaccard bundles (threshold = 0.5)

Arrays grouped by the seed-and-grow algorithm: every pair within a bundle has J >= 0.5. Seed order is by descending nest count.

16 bundles: sizes [3, 1, 1, 2, 1, 1, 1, 3, 5, 7, 2, 1, 1, 3, 1, 1], total 34 arrays.

### Bundle 1 (3 arrays) — J: min=0.50 avg=0.56 max=0.67

- Nests ∩: [1, 2, 4, 17]
- Nests ∪: [1, 2, 3, 4, 5, 6, 7, 9, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%vt` | 7 | [1, 2, 3, 4, 5, 6, 17] |
| `p_diag%vn_ie` | 5 | [1, 2, 4, 7, 17] |
| `z_kin_hor_e` | 5 | [1, 2, 4, 9, 17] |

### Bundle 2 (1 arrays) — J: —

- Nests ∩: [1, 2, 3, 4, 6]
- Nests ∪: [1, 2, 3, 4, 6]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_prog%vn` | 5 | [1, 2, 3, 4, 6] |

### Bundle 3 (1 arrays) — J: —

- Nests ∩: [8, 12, 13, 14, 15]
- Nests ∪: [8, 12, 13, 14, 15]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_con_c` | 5 | [8, 12, 13, 14, 15] |

### Bundle 4 (2 arrays) — J: min=0.50 avg=0.50 max=0.50

- Nests ∩: [1, 2]
- Nests ∪: [1, 2, 5, 7]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_vt_ie` | 4 | [1, 2, 5, 7] |
| `p_metrics%wgtfacq_e` | 2 | [1, 2] |

### Bundle 5 (1 arrays) — J: —

- Nests ∩: [7, 12, 15]
- Nests ∪: [7, 12, 15]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_prog%w` | 3 | [7, 12, 15] |

### Bundle 6 (1 arrays) — J: —

- Nests ∩: [9, 10, 16]
- Nests ∪: [9, 10, 16]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_int%e_bln_c_s` | 3 | [9, 10, 16] |

### Bundle 7 (1 arrays) — J: —

- Nests ∩: [4, 5]
- Nests ∪: [4, 5]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_metrics%wgtfac_e` | 2 | [4, 5] |

### Bundle 8 (3 arrays) — J: min=0.50 avg=0.67 max=1.00

- Nests ∩: [6]
- Nests ∪: [6, 10]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_concorr_me` | 2 | [6, 10] |
| `p_metrics%ddxn_z_full` | 1 | [6] |
| `p_metrics%ddxt_z_full` | 1 | [6] |

### Bundle 9 (5 arrays) — J: min=0.50 avg=0.80 max=1.00

- Nests ∩: [7]
- Nests ∪: [7, 16]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_v_grad_w` | 2 | [7, 16] |
| `p_patch%edges%inv_dual_edge_length` | 1 | [7] |
| `p_patch%edges%inv_primal_edge_length` | 1 | [7] |
| `p_patch%edges%tangent_orientation` | 1 | [7] |
| `z_w_v` | 1 | [7] |

### Bundle 10 (7 arrays) — J: min=0.50 avg=0.86 max=1.00

- Nests ∩: [17]
- Nests ∪: [9, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_ekinh` | 2 | [9, 17] |
| `p_diag%ddt_vn_apc_pc` | 1 | [17] |
| `p_int%c_lin_e` | 1 | [17] |
| `p_metrics%coeff_gradekin` | 1 | [17] |
| `p_metrics%ddqz_z_full_e` | 1 | [17] |
| `p_patch%edges%f_e` | 1 | [17] |
| `zeta` | 1 | [17] |

### Bundle 11 (2 arrays) — J: min=0.50 avg=0.50 max=0.50

- Nests ∩: [11]
- Nests ∪: [10, 11]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_concorr_mc` | 2 | [10, 11] |
| `p_metrics%wgtfac_c` | 1 | [11] |

### Bundle 12 (1 arrays) — J: —

- Nests ∩: [11, 13]
- Nests ∪: [11, 13]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%w_concorr_c` | 2 | [11, 13] |

### Bundle 13 (1 arrays) — J: —

- Nests ∩: [14, 17]
- Nests ∪: [14, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_con_c_full` | 2 | [14, 17] |

### Bundle 14 (3 arrays) — J: min=0.50 avg=0.67 max=1.00

- Nests ∩: [15]
- Nests ∪: [15, 16]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%ddt_w_adv_pc` | 2 | [15, 16] |
| `p_metrics%coeff1_dwdz` | 1 | [15] |
| `p_metrics%coeff2_dwdz` | 1 | [15] |

### Bundle 15 (1 arrays) — J: —

- Nests ∩: [2]
- Nests ∪: [2]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%vn_ie_ubc` | 1 | [2] |

### Bundle 16 (1 arrays) — J: —

- Nests ∩: [3]
- Nests ∪: [3]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_int%rbf_vec_coeff_e` | 1 | [3] |

## 9. Size-constrained bundles (valid sizes: 2, 4, 8, multiples of 8)

The unconstrained bundles from §8 are post-processed: invalid-sized bundles are merged (preferring highest average Jaccard) or padded (stealing the most-similar member from a donor) until every bundle reaches a valid size.

10 bundles: sizes [4, 2, 2, 4, 8, 4, 2, 2, 4, 2], total 34 arrays. All sizes valid.

### Bundle 1 (4 arrays) — J: min=0.43 avg=0.54 max=0.71

- Nests ∩: [1, 2, 4]
- Nests ∪: [1, 2, 3, 4, 5, 6, 7, 9, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_prog%vn` | 5 | [1, 2, 3, 4, 6] |
| `p_diag%vt` | 7 | [1, 2, 3, 4, 5, 6, 17] |
| `p_diag%vn_ie` | 5 | [1, 2, 4, 7, 17] |
| `z_kin_hor_e` | 5 | [1, 2, 4, 9, 17] |

### Bundle 2 (2 arrays) — J: min=0.33 avg=0.33 max=0.33

- Nests ∩: [12, 15]
- Nests ∪: [7, 8, 12, 13, 14, 15]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_con_c` | 5 | [8, 12, 13, 14, 15] |
| `p_prog%w` | 3 | [7, 12, 15] |

### Bundle 3 (2 arrays) — J: min=0.50 avg=0.50 max=0.50

- Nests ∩: [1, 2]
- Nests ∪: [1, 2, 5, 7]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_vt_ie` | 4 | [1, 2, 5, 7] |
| `p_metrics%wgtfacq_e` | 2 | [1, 2] |

### Bundle 4 (4 arrays) — J: min=0.00 avg=0.38 max=1.00

- Nests ∩: (none)
- Nests ∪: [6, 9, 10, 16]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_int%e_bln_c_s` | 3 | [9, 10, 16] |
| `z_w_concorr_me` | 2 | [6, 10] |
| `p_metrics%ddxn_z_full` | 1 | [6] |
| `p_metrics%ddxt_z_full` | 1 | [6] |

### Bundle 5 (8 arrays) — J: min=0.00 avg=0.64 max=1.00

- Nests ∩: (none)
- Nests ∪: [4, 5, 9, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_metrics%wgtfac_e` | 2 | [4, 5] |
| `z_ekinh` | 2 | [9, 17] |
| `p_diag%ddt_vn_apc_pc` | 1 | [17] |
| `p_int%c_lin_e` | 1 | [17] |
| `p_metrics%coeff_gradekin` | 1 | [17] |
| `p_metrics%ddqz_z_full_e` | 1 | [17] |
| `p_patch%edges%f_e` | 1 | [17] |
| `zeta` | 1 | [17] |

### Bundle 6 (4 arrays) — J: min=1.00 avg=1.00 max=1.00

- Nests ∩: [7]
- Nests ∪: [7]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_patch%edges%inv_dual_edge_length` | 1 | [7] |
| `p_patch%edges%inv_primal_edge_length` | 1 | [7] |
| `p_patch%edges%tangent_orientation` | 1 | [7] |
| `z_w_v` | 1 | [7] |

### Bundle 7 (2 arrays) — J: min=0.50 avg=0.50 max=0.50

- Nests ∩: [11]
- Nests ∪: [10, 11]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_concorr_mc` | 2 | [10, 11] |
| `p_metrics%wgtfac_c` | 1 | [11] |

### Bundle 8 (2 arrays) — J: min=0.00 avg=0.00 max=0.00

- Nests ∩: (none)
- Nests ∪: [11, 13, 14, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%w_concorr_c` | 2 | [11, 13] |
| `z_w_con_c_full` | 2 | [14, 17] |

### Bundle 9 (4 arrays) — J: min=0.00 avg=0.33 max=1.00

- Nests ∩: (none)
- Nests ∪: [2, 15, 16]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%vn_ie_ubc` | 1 | [2] |
| `p_diag%ddt_w_adv_pc` | 2 | [15, 16] |
| `p_metrics%coeff1_dwdz` | 1 | [15] |
| `p_metrics%coeff2_dwdz` | 1 | [15] |

### Bundle 10 (2 arrays) — J: min=0.00 avg=0.00 max=0.00

- Nests ∩: (none)
- Nests ∪: [3, 7, 16]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_int%rbf_vec_coeff_e` | 1 | [3] |
| `z_v_grad_w` | 2 | [7, 16] |

## 10. Bundle comparison summary

| | Unconstrained | Constrained |
|---|---:|---:|
| Bundles | 16 | 10 |
| Sizes | [3, 1, 1, 2, 1, 1, 1, 3, 5, 7, 2, 1, 1, 3, 1, 1] | [4, 2, 2, 4, 8, 4, 2, 2, 4, 2] |
| Singletons | 9 | 0 |
| Largest | 7 | 8 |
| Valid-sized | 2/16 | 10/10 |

## 11. Utilization-aware bundles (zero cache waste)

A bundle packed as AoS or AoSoA wastes cache-line bandwidth whenever a nest accesses some but not all members.  The **bundle utilization** measures the worst case:

$$\text{util}(B) = \min_{N:\, N \cap B \neq \emptyset} \frac{|N \cap B|}{|B|}$$

util${}=1.0$: every nest uses every field (zero waste). util${}=0.5$: some nest loads the record but ignores half.

### Utilization of Jaccard bundles (§9)

| Bundle | Size | Util | Worst nest | Accessed | Wasted |
|-------:|-----:|-----:|-----------:|---------:|-------:|
| 1 | 4 | 25% ⚠ | N5 | 1/4 | `p_prog%vn`, `p_diag%vn_ie`, `z_kin_hor_e` |
| 2 | 2 | 50% ⚠ | N7 | 1/2 | `z_w_con_c` |
| 3 | 2 | 50% ⚠ | N5 | 1/2 | `p_metrics%wgtfacq_e` |
| 4 | 4 | 25% ⚠ | N16 | 1/4 | `z_w_concorr_me`, `p_metrics%ddxn_z_full`, `p_metrics%ddxt_z_full` |
| 5 | 8 | 12% ⚠ | N4 | 1/8 | `z_ekinh`, `p_diag%ddt_vn_apc_pc`, `p_int%c_lin_e`, `p_metrics%coeff_gradekin`, `p_metrics%ddqz_z_full_e`, `p_patch%edges%f_e`, `zeta` |
| 6 | 4 | 100% | — | 0/4 | — |
| 7 | 2 | 50% ⚠ | N10 | 1/2 | `p_metrics%wgtfac_c` |
| 8 | 2 | 50% ⚠ | N17 | 1/2 | `p_diag%w_concorr_c` |
| 9 | 4 | 25% ⚠ | N16 | 1/4 | `p_diag%vn_ie_ubc`, `p_metrics%coeff1_dwdz`, `p_metrics%coeff2_dwdz` |
| 10 | 2 | 50% ⚠ | N16 | 1/2 | `p_int%rbf_vec_coeff_e` |

### Strict co-access bundles (util = 1.0)

Arrays with identical nest-membership sets are grouped, then groups are merged only if the merged bundle maintains util = 1.0 (every nest that touches the bundle uses every field).

24 bundles: sizes [6, 4, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], total 34 arrays.

### Strict bundle 1 (6 arrays) — util=100%

- Nests ∩: [17]
- Nests ∪: [17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%ddt_vn_apc_pc` | 1 | [17] |
| `p_int%c_lin_e` | 1 | [17] |
| `p_metrics%coeff_gradekin` | 1 | [17] |
| `p_metrics%ddqz_z_full_e` | 1 | [17] |
| `p_patch%edges%f_e` | 1 | [17] |
| `zeta` | 1 | [17] |

### Strict bundle 2 (4 arrays) — util=100%

- Nests ∩: [7]
- Nests ∪: [7]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_patch%edges%inv_dual_edge_length` | 1 | [7] |
| `p_patch%edges%inv_primal_edge_length` | 1 | [7] |
| `p_patch%edges%tangent_orientation` | 1 | [7] |
| `z_w_v` | 1 | [7] |

### Strict bundle 3 (2 arrays) — util=100%

- Nests ∩: [15]
- Nests ∪: [15]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_metrics%coeff1_dwdz` | 1 | [15] |
| `p_metrics%coeff2_dwdz` | 1 | [15] |

### Strict bundle 4 (2 arrays) — util=100%

- Nests ∩: [6]
- Nests ∪: [6]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_metrics%ddxn_z_full` | 1 | [6] |
| `p_metrics%ddxt_z_full` | 1 | [6] |

### Strict bundle 5 (1 arrays) — util=100%

- Nests ∩: [15, 16]
- Nests ∪: [15, 16]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%ddt_w_adv_pc` | 2 | [15, 16] |

### Strict bundle 6 (1 arrays) — util=100%

- Nests ∩: [1, 2, 4, 7, 17]
- Nests ∪: [1, 2, 4, 7, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%vn_ie` | 5 | [1, 2, 4, 7, 17] |

### Strict bundle 7 (1 arrays) — util=100%

- Nests ∩: [2]
- Nests ∪: [2]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%vn_ie_ubc` | 1 | [2] |

### Strict bundle 8 (1 arrays) — util=100%

- Nests ∩: [1, 2, 3, 4, 5, 6, 17]
- Nests ∪: [1, 2, 3, 4, 5, 6, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%vt` | 7 | [1, 2, 3, 4, 5, 6, 17] |

### Strict bundle 9 (1 arrays) — util=100%

- Nests ∩: [11, 13]
- Nests ∪: [11, 13]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_diag%w_concorr_c` | 2 | [11, 13] |

### Strict bundle 10 (1 arrays) — util=100%

- Nests ∩: [9, 10, 16]
- Nests ∪: [9, 10, 16]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_int%e_bln_c_s` | 3 | [9, 10, 16] |

### Strict bundle 11 (1 arrays) — util=100%

- Nests ∩: [3]
- Nests ∪: [3]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_int%rbf_vec_coeff_e` | 1 | [3] |

### Strict bundle 12 (1 arrays) — util=100%

- Nests ∩: [11]
- Nests ∪: [11]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_metrics%wgtfac_c` | 1 | [11] |

### Strict bundle 13 (1 arrays) — util=100%

- Nests ∩: [4, 5]
- Nests ∪: [4, 5]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_metrics%wgtfac_e` | 2 | [4, 5] |

### Strict bundle 14 (1 arrays) — util=100%

- Nests ∩: [1, 2]
- Nests ∪: [1, 2]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_metrics%wgtfacq_e` | 2 | [1, 2] |

### Strict bundle 15 (1 arrays) — util=100%

- Nests ∩: [1, 2, 3, 4, 6]
- Nests ∪: [1, 2, 3, 4, 6]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_prog%vn` | 5 | [1, 2, 3, 4, 6] |

### Strict bundle 16 (1 arrays) — util=100%

- Nests ∩: [7, 12, 15]
- Nests ∪: [7, 12, 15]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `p_prog%w` | 3 | [7, 12, 15] |

### Strict bundle 17 (1 arrays) — util=100%

- Nests ∩: [9, 17]
- Nests ∪: [9, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_ekinh` | 2 | [9, 17] |

### Strict bundle 18 (1 arrays) — util=100%

- Nests ∩: [1, 2, 4, 9, 17]
- Nests ∪: [1, 2, 4, 9, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_kin_hor_e` | 5 | [1, 2, 4, 9, 17] |

### Strict bundle 19 (1 arrays) — util=100%

- Nests ∩: [7, 16]
- Nests ∪: [7, 16]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_v_grad_w` | 2 | [7, 16] |

### Strict bundle 20 (1 arrays) — util=100%

- Nests ∩: [1, 2, 5, 7]
- Nests ∪: [1, 2, 5, 7]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_vt_ie` | 4 | [1, 2, 5, 7] |

### Strict bundle 21 (1 arrays) — util=100%

- Nests ∩: [8, 12, 13, 14, 15]
- Nests ∪: [8, 12, 13, 14, 15]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_con_c` | 5 | [8, 12, 13, 14, 15] |

### Strict bundle 22 (1 arrays) — util=100%

- Nests ∩: [14, 17]
- Nests ∪: [14, 17]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_con_c_full` | 2 | [14, 17] |

### Strict bundle 23 (1 arrays) — util=100%

- Nests ∩: [10, 11]
- Nests ∪: [10, 11]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_concorr_mc` | 2 | [10, 11] |

### Strict bundle 24 (1 arrays) — util=100%

- Nests ∩: [6, 10]
- Nests ∪: [6, 10]

| Array | #Nests | Nests |
|-------|:------:|-------|
| `z_w_concorr_me` | 2 | [6, 10] |

### Strict bundles with size constraint

13 bundles: sizes [8, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

| Bundle | Size | Util | Arrays |
|-------:|-----:|-----:|--------|
| 1 | 8 | 12% ⚠ | `p_diag%ddt_vn_apc_pc`, `p_int%c_lin_e`, `p_metrics%coeff_gradekin`, `p_metrics%ddqz_z_full_e`, `p_patch%edges%f_e`, `zeta`, `p_int%e_bln_c_s`, `z_ekinh` |
| 2 | 4 | 100% | `p_patch%edges%inv_dual_edge_length`, `p_patch%edges%inv_primal_edge_length`, `p_patch%edges%tangent_orientation`, `z_w_v` |
| 3 | 2 | 100% | `p_metrics%coeff1_dwdz`, `p_metrics%coeff2_dwdz` |
| 4 | 2 | 100% | `p_metrics%ddxn_z_full`, `p_metrics%ddxt_z_full` |
| 5 | 2 | 50% ⚠ | `p_diag%ddt_w_adv_pc`, `z_v_grad_w` |
| 6 | 2 | 50% ⚠ | `p_diag%vn_ie`, `z_kin_hor_e` |
| 7 | 2 | 50% ⚠ | `p_diag%vn_ie_ubc`, `p_metrics%wgtfacq_e` |
| 8 | 2 | 50% ⚠ | `p_diag%vt`, `p_prog%vn` |
| 9 | 2 | 50% ⚠ | `p_diag%w_concorr_c`, `p_metrics%wgtfac_c` |
| 10 | 2 | 50% ⚠ | `p_int%rbf_vec_coeff_e`, `p_metrics%wgtfac_e` |
| 11 | 2 | 50% ⚠ | `p_prog%w`, `z_w_con_c` |
| 12 | 2 | 50% ⚠ | `z_vt_ie`, `z_w_con_c_full` |
| 13 | 2 | 50% ⚠ | `z_w_concorr_mc`, `z_w_concorr_me` |

### Utilization comparison

| Method | Bundles | Min util | Avg util | Bundles with util<1 |
|--------|--------:|---------:|---------:|--------------------:|
| Jaccard §9 | 10 | 12% | 44% | 9 |
| Strict | 24 | 100% | 100% | 0 |
| Strict+sized | 13 | 12% | 59% | 10 |

