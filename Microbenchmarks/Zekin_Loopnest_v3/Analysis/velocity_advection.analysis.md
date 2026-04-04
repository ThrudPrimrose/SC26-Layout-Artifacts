# Loop and Array Access Analysis: `velocity_tendencies`

| Symbol | Meaning |
|--------|---------|
| `h` | horizontal (nproma) dimension — loop vars `jc`, `je` (treated as equivalent) |
| `v` | vertical (nlev) dimension — loop var `jk` |
| `b` | block dimension — loop var `jb` |
| `S` | **structured** access — subscript is a direct loop variable |
| `U` | **unstructured** access — subscript wraps an indirection array (`icidx`, `ividx`, ...) |
| `C` | **constant** access — subscript is loop-invariant |

## Array Groups

Arrays in the same group have identical access patterns across all loops.

### FULLY STRUCTURED (8 arrays)

Patterns observed:

- `h:S  v:S  b:S`

Arrays:

- `p_diag%w_concorr_c`
- `p_metrics%coeff1_dwdz`
- `p_metrics%coeff2_dwdz`
- `p_metrics%ddqz_z_full_e`
- `p_metrics%ddxn_z_full`
- `p_metrics%ddxt_z_full`
- `p_metrics%wgtfac_c`
- `p_metrics%wgtfac_e`

### FULLY STRUCTURED (5 arrays)

Patterns observed:

- `h:S  -:C  b:S`

Arrays:

- `p_diag%vn_ie_ubc`
- `p_int%c_lin_e`
- `p_int%e_bln_c_s`
- `p_metrics%coeff_gradekin`
- `p_metrics%wgtfacq_e`

### HAS UNSTRUCTURED (5 arrays)

Patterns observed:

- `h:S  v:S  b:S`
- `h:U  v:S  h:U`

Arrays:

- `p_prog%w`
- `z_ekinh`
- `z_v_grad_w`
- `z_w_con_c_full`
- `z_w_concorr_me`

### FULLY STRUCTURED (4 arrays)

Patterns observed:

- `h:S  b:S`

Arrays:

- `p_patch%edges%f_e`
- `p_patch%edges%inv_dual_edge_length`
- `p_patch%edges%inv_primal_edge_length`
- `p_patch%edges%tangent_orientation`

### FULLY STRUCTURED (3 arrays)

Patterns observed:

- `h:S  -:C  b:S`
- `h:S  v:S  b:S`

Arrays:

- `p_diag%vn_ie`
- `p_diag%vt`
- `z_vt_ie`

### HAS UNSTRUCTURED (2 arrays)

Patterns observed:

- `h:S  -:C  b:S`
- `h:S  v:S  b:S`
- `h:U  v:S  h:U`

Arrays:

- `p_prog%vn`
- `z_kin_hor_e`

### HAS UNSTRUCTURED (2 arrays)

Patterns observed:

- `h:U  v:S  h:U`

Arrays:

- `z_w_v`
- `zeta`

### FULLY STRUCTURED (2 arrays)

Patterns observed:

- `h:S  v:S  b:S  -:C`

Arrays:

- `p_diag%ddt_vn_apc_pc`
- `p_diag%ddt_w_adv_pc`

### FULLY STRUCTURED (1 arrays)

Patterns observed:

- `-:C  h:S  b:S`

Arrays:

- `p_int%rbf_vec_coeff_e`

### FULLY STRUCTURED (1 arrays)

Patterns observed:

- `h:S  -:C`
- `h:S  v:S`

Arrays:

- `z_w_con_c`

### FULLY STRUCTURED (1 arrays)

Patterns observed:

- `h:S  v:S`

Arrays:

- `z_w_concorr_mc`

## Per-Array Access Summary

> `S+U` marks arrays accessed both structured and unstructured (layout tradeoff)

| Array | Role:S/U Patterns | Conflict |
|-------|-------------------|----------|
| `p_diag%ddt_vn_apc_pc` | `h:S  v:S  b:S  -:C` |  |
| `p_diag%ddt_w_adv_pc` | `h:S  v:S  b:S  -:C` |  |
| `p_diag%vn_ie` | `h:S  -:C  b:S` \| `h:S  v:S  b:S` |  |
| `p_diag%vn_ie_ubc` | `h:S  -:C  b:S` |  |
| `p_diag%vt` | `h:S  -:C  b:S` \| `h:S  v:S  b:S` |  |
| `p_diag%w_concorr_c` | `h:S  v:S  b:S` |  |
| `p_int%c_lin_e` | `h:S  -:C  b:S` |  |
| `p_int%e_bln_c_s` | `h:S  -:C  b:S` |  |
| `p_int%rbf_vec_coeff_e` | `-:C  h:S  b:S` |  |
| `p_metrics%coeff1_dwdz` | `h:S  v:S  b:S` |  |
| `p_metrics%coeff2_dwdz` | `h:S  v:S  b:S` |  |
| `p_metrics%coeff_gradekin` | `h:S  -:C  b:S` |  |
| `p_metrics%ddqz_z_full_e` | `h:S  v:S  b:S` |  |
| `p_metrics%ddxn_z_full` | `h:S  v:S  b:S` |  |
| `p_metrics%ddxt_z_full` | `h:S  v:S  b:S` |  |
| `p_metrics%wgtfac_c` | `h:S  v:S  b:S` |  |
| `p_metrics%wgtfac_e` | `h:S  v:S  b:S` |  |
| `p_metrics%wgtfacq_e` | `h:S  -:C  b:S` |  |
| `p_patch%edges%f_e` | `h:S  b:S` |  |
| `p_patch%edges%inv_dual_edge_length` | `h:S  b:S` |  |
| `p_patch%edges%inv_primal_edge_length` | `h:S  b:S` |  |
| `p_patch%edges%tangent_orientation` | `h:S  b:S` |  |
| `p_prog%vn` | `h:S  -:C  b:S` \| `h:S  v:S  b:S` \| `h:U  v:S  h:U` |  |
| `p_prog%w` | `h:S  v:S  b:S` \| `h:U  v:S  h:U` |  |
| `z_ekinh` | `h:S  v:S  b:S` \| `h:U  v:S  h:U` |  |
| `z_kin_hor_e` | `h:S  -:C  b:S` \| `h:S  v:S  b:S` \| `h:U  v:S  h:U` |  |
| `z_v_grad_w` | `h:S  v:S  b:S` \| `h:U  v:S  h:U` |  |
| `z_vt_ie` | `h:S  -:C  b:S` \| `h:S  v:S  b:S` |  |
| `z_w_con_c` | `h:S  -:C` \| `h:S  v:S` |  |
| `z_w_con_c_full` | `h:S  v:S  b:S` \| `h:U  v:S  h:U` |  |
| `z_w_concorr_mc` | `h:S  v:S` |  |
| `z_w_concorr_me` | `h:S  v:S  b:S` \| `h:U  v:S  h:U` |  |
| `z_w_v` | `h:U  v:S  h:U` |  |
| `zeta` | `h:U  v:S  h:U` |  |

## Loop Shape Counts

| Count | Shape | Ranges | Behavior |
|------:|-------|--------|----------|
| 6 | `v.h` | full_vert, full_horiz |  |
| 6 | `v.h` | partial_vert, full_horiz |  |
| 2 | `b.h` | full_block, full_horiz |  |
| 2 | `v.h` | partial_vert, full_horiz | accumulate |
| 1 | `b.h` | full_block, full_horiz | accumulate |

## Loop Patterns (collapsed)

> **Collapsed**: if ANY array in the nest is `U` for a role, the whole nest is `U` for that role.

| Count | Shape | Ranges | Behavior | Collapsed S/U |
|------:|-------|--------|----------|---------------|
| 5 | `v.h` | partial_vert, full_horiz |  | `h:S  v:S  b:S` |
| 4 | `v.h` | full_vert, full_horiz |  | `h:U  v:S  b:S` |
| 2 | `v.h` | full_vert, full_horiz |  | `h:S  v:S  b:S` |
| 1 | `b.h` | full_block, full_horiz |  | `h:S  b:S` |
| 1 | `b.h` | full_block, full_horiz | accumulate | `h:S  b:S` |
| 1 | `b.h` | full_block, full_horiz |  | `h:S` |
| 1 | `v.h` | partial_vert, full_horiz |  | `h:U  v:S  b:S` |
| 1 | `v.h` | partial_vert, full_horiz | accumulate | `h:S  v:S  b:S` |
| 1 | `v.h` | partial_vert, full_horiz | accumulate | `h:U  v:S  b:S` |

## Detailed Loop Nests

17 nests found.

### Nest 1: `b.h` (full_block, full_horiz)

Collapsed: `h:S  b:S`

```fortran
! inside DO jb = i_startblk:i_endblk
DO je=i_startidx,i_endidx
  p_diag%vn_ie(je, 1, jb) = p_prog%vn(je, 1, jb)
  z_vt_ie(je, 1, jb) = p_diag%vt(je, 1, jb)
  z_kin_hor_e(je, 1, jb) = 0.5_wp*(p_prog%vn(je, 1, jb)**2 + p_diag%vt(je, 1, jb)**2)
  p_diag%vn_ie(je, nlevp1, jb) = p_metrics%wgtfacq_e(je, 1, jb)*p_prog%vn(je, nlev, jb) + p_metrics%wgtfacq_e(je, 2, jb) &
  & *p_prog%vn(je, nlev - 1, jb) + p_metrics%wgtfacq_e(je, 3, jb)*p_prog%vn(je, nlev - 2, jb)
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `p_diag%vn_ie` | `h:S  -:C  b:S` |
| `p_prog%vn` | `h:S  -:C  b:S` |
| `z_vt_ie` | `h:S  -:C  b:S` |
| `p_diag%vt` | `h:S  -:C  b:S` |
| `z_kin_hor_e` | `h:S  -:C  b:S` |
| `p_metrics%wgtfacq_e` | `h:S  -:C  b:S` |

</details>

### Nest 2: `b.h` (full_block, full_horiz) `accumulate`

Collapsed: `h:S  b:S`

```fortran
! inside DO jb = i_startblk:i_endblk
DO je=i_startidx,i_endidx
  p_diag%vn_ie(je, 1, jb) = p_diag%vn_ie_ubc(je, 1, jb) + dt_linintp_ubc*p_diag%vn_ie_ubc(je, 2, jb)
  z_vt_ie(je, 1, jb) = p_diag%vt(je, 1, jb)
  z_kin_hor_e(je, 1, jb) = 0.5_wp*(p_prog%vn(je, 1, jb)**2 + p_diag%vt(je, 1, jb)**2)
  p_diag%vn_ie(je, nlevp1, jb) = p_metrics%wgtfacq_e(je, 1, jb)*p_prog%vn(je, nlev, jb) + p_metrics%wgtfacq_e(je, 2, jb) &
  & *p_prog%vn(je, nlev - 1, jb) + p_metrics%wgtfacq_e(je, 3, jb)*p_prog%vn(je, nlev - 2, jb)
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `p_diag%vn_ie` | `h:S  -:C  b:S` |
| `p_diag%vn_ie_ubc` | `h:S  -:C  b:S` |
| `z_vt_ie` | `h:S  -:C  b:S` |
| `p_diag%vt` | `h:S  -:C  b:S` |
| `z_kin_hor_e` | `h:S  -:C  b:S` |
| `p_prog%vn` | `h:S  -:C  b:S` |
| `p_metrics%wgtfacq_e` | `h:S  -:C  b:S` |

</details>

### Nest 3: `v.h` (full_vert, full_horiz)

Collapsed: `h:U  v:S  b:S`

```fortran
DO jk=1,nlev
  DO je=i_startidx,i_endidx
    p_diag%vt(je, jk, jb) = p_int%rbf_vec_coeff_e(1, je, jb)*p_prog%vn(iqidx(je, jb, 1), jk, iqblk(je, jb, 1)) +  &
    & p_int%rbf_vec_coeff_e(2, je, jb)*p_prog%vn(iqidx(je, jb, 2), jk, iqblk(je, jb, 2)) + p_int%rbf_vec_coeff_e(3, je, jb) &
    & *p_prog%vn(iqidx(je, jb, 3), jk, iqblk(je, jb, 3)) + p_int%rbf_vec_coeff_e(4, je, jb)*p_prog%vn(iqidx(je, jb, 4), jk,  &
    & iqblk(je, jb, 4))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `p_diag%vt` | `h:S  v:S  b:S` |
| `p_int%rbf_vec_coeff_e` | `-:C  h:S  b:S` |
| `p_prog%vn` | `h:U  v:S  h:U` |

</details>

### Nest 4: `v.h` (partial_vert, full_horiz)

Collapsed: `h:S  v:S  b:S`

```fortran
DO jk=2,nlev
  DO je=i_startidx,i_endidx
    p_diag%vn_ie(je, jk, jb) =  &
    & p_metrics%wgtfac_e(je, jk, jb)*p_prog%vn(je, jk, jb) + (1._wp - p_metrics%wgtfac_e(je, jk, jb))*p_prog%vn(je, jk - 1, jb)
    z_kin_hor_e(je, jk, jb) = 0.5_wp*(p_prog%vn(je, jk, jb)**2 + p_diag%vt(je, jk, jb)**2)
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `p_diag%vn_ie` | `h:S  v:S  b:S` |
| `p_metrics%wgtfac_e` | `h:S  v:S  b:S` |
| `p_prog%vn` | `h:S  v:S  b:S` |
| `z_kin_hor_e` | `h:S  v:S  b:S` |
| `p_diag%vt` | `h:S  v:S  b:S` |

</details>

### Nest 5: `v.h` (partial_vert, full_horiz)

Collapsed: `h:S  v:S  b:S`

```fortran
DO jk=2,nlev
  DO je=i_startidx,i_endidx
    z_vt_ie(je, jk, jb) =  &
    & p_metrics%wgtfac_e(je, jk, jb)*p_diag%vt(je, jk, jb) + (1._wp - p_metrics%wgtfac_e(je, jk, jb))*p_diag%vt(je, jk - 1, jb)
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `z_vt_ie` | `h:S  v:S  b:S` |
| `p_metrics%wgtfac_e` | `h:S  v:S  b:S` |
| `p_diag%vt` | `h:S  v:S  b:S` |

</details>

### Nest 6: `v.h` (partial_vert, full_horiz)

Collapsed: `h:S  v:S  b:S`

```fortran
DO jk=nflatlev_jg,nlev
  DO je=i_startidx,i_endidx
    z_w_concorr_me(je, jk, jb) =  &
    & p_prog%vn(je, jk, jb)*p_metrics%ddxn_z_full(je, jk, jb) + p_diag%vt(je, jk, jb)*p_metrics%ddxt_z_full(je, jk, jb)
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `z_w_concorr_me` | `h:S  v:S  b:S` |
| `p_prog%vn` | `h:S  v:S  b:S` |
| `p_metrics%ddxn_z_full` | `h:S  v:S  b:S` |
| `p_diag%vt` | `h:S  v:S  b:S` |
| `p_metrics%ddxt_z_full` | `h:S  v:S  b:S` |

</details>

### Nest 7: `v.h` (full_vert, full_horiz)

Collapsed: `h:U  v:S  b:S`

```fortran
DO jk=1,nlev
  DO je=i_startidx,i_endidx
    z_v_grad_w(je, jk, jb) = p_diag%vn_ie(je, jk, jb)*p_patch%edges%inv_dual_edge_length(je, jb)*(p_prog%w(icidx(je, jb, 1), jk,  &
    & icblk(je, jb, 1)) - p_prog%w(icidx(je, jb, 2), jk, icblk(je, jb, 2))) + z_vt_ie(je, jk, jb) &
    & *p_patch%edges%inv_primal_edge_length(je, jb)*p_patch%edges%tangent_orientation(je, jb)*(z_w_v(ividx(je, jb, 1), jk,  &
    & ivblk(je, jb, 1)) - z_w_v(ividx(je, jb, 2), jk, ivblk(je, jb, 2)))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `z_v_grad_w` | `h:S  v:S  b:S` |
| `p_diag%vn_ie` | `h:S  v:S  b:S` |
| `p_patch%edges%inv_dual_edge_length` | `h:S  b:S` |
| `p_prog%w` | `h:U  v:S  h:U` |
| `z_vt_ie` | `h:S  v:S  b:S` |
| `p_patch%edges%inv_primal_edge_length` | `h:S  b:S` |
| `p_patch%edges%tangent_orientation` | `h:S  b:S` |
| `z_w_v` | `h:U  v:S  h:U` |

</details>

### Nest 8: `b.h` (full_block, full_horiz)

Collapsed: `h:S`

```fortran
! inside DO jb = i_startblk:i_endblk
DO jc=i_startidx,i_endidx
  z_w_con_c(jc, nlevp1) = 0.0_vp
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `z_w_con_c` | `h:S  -:C` |

</details>

### Nest 9: `v.h` (full_vert, full_horiz)

Collapsed: `h:U  v:S  b:S`

```fortran
DO jk=1,nlev
  DO jc=i_startidx,i_endidx
    z_ekinh(jc, jk, jb) = p_int%e_bln_c_s(jc, 1, jb)*z_kin_hor_e(ieidx(jc, jb, 1), jk, ieblk(jc, jb, 1)) + p_int%e_bln_c_s(jc, 2, &
    &  jb)*z_kin_hor_e(ieidx(jc, jb, 2), jk, ieblk(jc, jb, 2)) + p_int%e_bln_c_s(jc, 3, jb)*z_kin_hor_e(ieidx(jc, jb, 3), jk,  &
    & ieblk(jc, jb, 3))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `z_ekinh` | `h:S  v:S  b:S` |
| `p_int%e_bln_c_s` | `h:S  -:C  b:S` |
| `z_kin_hor_e` | `h:U  v:S  h:U` |

</details>

### Nest 10: `v.h` (partial_vert, full_horiz)

Collapsed: `h:U  v:S  b:S`

```fortran
DO jk=nflatlev_jg,nlev
  DO jc=i_startidx,i_endidx
    z_w_concorr_mc(jc, jk) = p_int%e_bln_c_s(jc, 1, jb)*z_w_concorr_me(ieidx(jc, jb, 1), jk, ieblk(jc, jb, 1)) +  &
    & p_int%e_bln_c_s(jc, 2, jb)*z_w_concorr_me(ieidx(jc, jb, 2), jk, ieblk(jc, jb, 2)) + p_int%e_bln_c_s(jc, 3, jb) &
    & *z_w_concorr_me(ieidx(jc, jb, 3), jk, ieblk(jc, jb, 3))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `z_w_concorr_mc` | `h:S  v:S` |
| `p_int%e_bln_c_s` | `h:S  -:C  b:S` |
| `z_w_concorr_me` | `h:U  v:S  h:U` |

</details>

### Nest 11: `v.h` (partial_vert, full_horiz)

Collapsed: `h:S  v:S  b:S`

```fortran
DO jk=nflatlev_jg + 1,nlev
  DO jc=i_startidx,i_endidx
    p_diag%w_concorr_c(jc, jk, jb) =  &
    & p_metrics%wgtfac_c(jc, jk, jb)*z_w_concorr_mc(jc, jk) + (1._vp - p_metrics%wgtfac_c(jc, jk, jb))*z_w_concorr_mc(jc, jk - 1)
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `p_diag%w_concorr_c` | `h:S  v:S  b:S` |
| `p_metrics%wgtfac_c` | `h:S  v:S  b:S` |
| `z_w_concorr_mc` | `h:S  v:S` |

</details>

### Nest 12: `v.h` (full_vert, full_horiz)

Collapsed: `h:S  v:S  b:S`

```fortran
DO jk=1,nlev
  DO jc=i_startidx,i_endidx
    z_w_con_c(jc, jk) = p_prog%w(jc, jk, jb)
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `z_w_con_c` | `h:S  v:S` |
| `p_prog%w` | `h:S  v:S  b:S` |

</details>

### Nest 13: `v.h` (partial_vert, full_horiz) `accumulate`

Collapsed: `h:S  v:S  b:S`

```fortran
DO jk=nlev,nflatlev_jg + 1,-1
  DO jc=i_startidx,i_endidx
    z_w_con_c(jc, jk) = z_w_con_c(jc, jk) - p_diag%w_concorr_c(jc, jk, jb)
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `z_w_con_c` | `h:S  v:S` |
| `p_diag%w_concorr_c` | `h:S  v:S  b:S` |

</details>

### Nest 14: `v.h` (full_vert, full_horiz)

Collapsed: `h:S  v:S  b:S`

```fortran
DO jk=1,nlev
  DO jc=i_startidx,i_endidx
    z_w_con_c_full(jc, jk, jb) = 0.5_vp*(z_w_con_c(jc, jk) + z_w_con_c(jc, jk + 1))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `z_w_con_c_full` | `h:S  v:S  b:S` |
| `z_w_con_c` | `h:S  v:S` |

</details>

### Nest 15: `v.h` (partial_vert, full_horiz)

Collapsed: `h:S  v:S  b:S`

```fortran
DO jk=2,nlev
  DO jc=i_startidx_2,i_endidx_2
    p_diag%ddt_w_adv_pc(jc, jk, jb, ntnd) = -z_w_con_c(jc, jk)*(p_prog%w(jc, jk - 1, jb)*p_metrics%coeff1_dwdz(jc, jk, jb) -  &
    & p_prog%w(jc, jk + 1, jb)*p_metrics%coeff2_dwdz(jc, jk, jb) + p_prog%w(jc, jk, jb)*(p_metrics%coeff2_dwdz(jc, jk, jb) -  &
    & p_metrics%coeff1_dwdz(jc, jk, jb)))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `p_diag%ddt_w_adv_pc` | `h:S  v:S  b:S  -:C` |
| `z_w_con_c` | `h:S  v:S` |
| `p_prog%w` | `h:S  v:S  b:S` |
| `p_metrics%coeff1_dwdz` | `h:S  v:S  b:S` |
| `p_metrics%coeff2_dwdz` | `h:S  v:S  b:S` |

</details>

### Nest 16: `v.h` (partial_vert, full_horiz) `accumulate`

Collapsed: `h:U  v:S  b:S`

```fortran
DO jk=2,nlev
  DO jc=i_startidx_2,i_endidx_2
    p_diag%ddt_w_adv_pc(jc, jk, jb, ntnd) = p_diag%ddt_w_adv_pc(jc, jk, jb, ntnd) + p_int%e_bln_c_s(jc, 1, jb) &
    & *z_v_grad_w(ieidx(jc, jb, 1), jk, ieblk(jc, jb, 1)) + p_int%e_bln_c_s(jc, 2, jb)*z_v_grad_w(ieidx(jc, jb, 2), jk, ieblk(jc, &
    &  jb, 2)) + p_int%e_bln_c_s(jc, 3, jb)*z_v_grad_w(ieidx(jc, jb, 3), jk, ieblk(jc, jb, 3))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `p_diag%ddt_w_adv_pc` | `h:S  v:S  b:S  -:C` |
| `p_int%e_bln_c_s` | `h:S  -:C  b:S` |
| `z_v_grad_w` | `h:U  v:S  h:U` |

</details>

### Nest 17: `v.h` (full_vert, full_horiz)

Collapsed: `h:U  v:S  b:S`

```fortran
DO jk=1,nlev
  DO je=i_startidx,i_endidx
    p_diag%ddt_vn_apc_pc(je, jk, jb, ntnd) = -(z_kin_hor_e(je, jk, jb)*(p_metrics%coeff_gradekin(je, 1, jb) -  &
    & p_metrics%coeff_gradekin(je, 2, jb)) + p_metrics%coeff_gradekin(je, 2, jb)*z_ekinh(icidx(je, jb, 2), jk, icblk(je, jb, 2))  &
    & - p_metrics%coeff_gradekin(je, 1, jb)*z_ekinh(icidx(je, jb, 1), jk, icblk(je, jb, 1)) + p_diag%vt(je, jk, jb) &
    & *(p_patch%edges%f_e(je, jb) + 0.5_vp*(zeta(ividx(je, jb, 1), jk, ivblk(je, jb, 1)) + zeta(ividx(je, jb, 2), jk, ivblk(je,  &
    & jb, 2)))) + (p_int%c_lin_e(je, 1, jb)*z_w_con_c_full(icidx(je, jb, 1), jk, icblk(je, jb, 1)) + p_int%c_lin_e(je, 2, jb) &
    & *z_w_con_c_full(icidx(je, jb, 2), jk, icblk(je, jb, 2)))*(p_diag%vn_ie(je, jk, jb) - p_diag%vn_ie(je, jk + 1, jb)) /  &
    & p_metrics%ddqz_z_full_e(je, jk, jb))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U |
|-------|----------|
| `p_diag%ddt_vn_apc_pc` | `h:S  v:S  b:S  -:C` |
| `z_kin_hor_e` | `h:S  v:S  b:S` |
| `p_metrics%coeff_gradekin` | `h:S  -:C  b:S` |
| `z_ekinh` | `h:U  v:S  h:U` |
| `p_diag%vt` | `h:S  v:S  b:S` |
| `p_patch%edges%f_e` | `h:S  b:S` |
| `zeta` | `h:U  v:S  h:U` |
| `p_int%c_lin_e` | `h:S  -:C  b:S` |
| `z_w_con_c_full` | `h:U  v:S  h:U` |
| `p_diag%vn_ie` | `h:S  v:S  b:S` |
| `p_metrics%ddqz_z_full_e` | `h:S  v:S  b:S` |

</details>


---

# Single-Block View (`nblks=1`)

Block dimension stripped: `b.h` becomes `h`, `b:S`/`b:U` removed from all signatures.

## Loop Patterns (single-block, collapsed)

| Count | Shape | Ranges | Behavior | Collapsed S/U |
|------:|-------|--------|----------|---------------|
| 5 | `v.h` | partial_vert, full_horiz |  | `h:S  v:S` |
| 4 | `v.h` | full_vert, full_horiz |  | `h:U  v:S` |
| 2 | `h` | full_horiz |  | `h:S` |
| 2 | `v.h` | full_vert, full_horiz |  | `h:S  v:S` |
| 1 | `h` | full_horiz | accumulate | `h:S` |
| 1 | `v.h` | partial_vert, full_horiz |  | `h:U  v:S` |
| 1 | `v.h` | partial_vert, full_horiz | accumulate | `h:S  v:S` |
| 1 | `v.h` | partial_vert, full_horiz | accumulate | `h:U  v:S` |

## Example Nest per Pattern Type

One representative nest from each pattern (the one touching the most arrays).

### `v.h` partial_vert, full_horiz — collapsed `h:S  v:S` (5x)

5 unique arrays.

```fortran
DO jk=2,nlev
  DO je=i_startidx,i_endidx
    p_diag%vn_ie(je, jk, jb) =  &
    & p_metrics%wgtfac_e(je, jk, jb)*p_prog%vn(je, jk, jb) + (1._wp - p_metrics%wgtfac_e(je, jk, jb))*p_prog%vn(je, jk - 1, jb)
    z_kin_hor_e(je, jk, jb) = 0.5_wp*(p_prog%vn(je, jk, jb)**2 + p_diag%vt(je, jk, jb)**2)
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U (no block) |
|-------|---------------------|
| `p_diag%vn_ie` | `h:S  v:S` |
| `p_metrics%wgtfac_e` | `h:S  v:S` |
| `p_prog%vn` | `h:S  v:S` |
| `z_kin_hor_e` | `h:S  v:S` |
| `p_diag%vt` | `h:S  v:S` |

</details>

### `v.h` full_vert, full_horiz — collapsed `h:U  v:S` (4x)

11 unique arrays.

```fortran
DO jk=1,nlev
  DO je=i_startidx,i_endidx
    p_diag%ddt_vn_apc_pc(je, jk, jb, ntnd) = -(z_kin_hor_e(je, jk, jb)*(p_metrics%coeff_gradekin(je, 1, jb) -  &
    & p_metrics%coeff_gradekin(je, 2, jb)) + p_metrics%coeff_gradekin(je, 2, jb)*z_ekinh(icidx(je, jb, 2), jk, icblk(je, jb, 2))  &
    & - p_metrics%coeff_gradekin(je, 1, jb)*z_ekinh(icidx(je, jb, 1), jk, icblk(je, jb, 1)) + p_diag%vt(je, jk, jb) &
    & *(p_patch%edges%f_e(je, jb) + 0.5_vp*(zeta(ividx(je, jb, 1), jk, ivblk(je, jb, 1)) + zeta(ividx(je, jb, 2), jk, ivblk(je,  &
    & jb, 2)))) + (p_int%c_lin_e(je, 1, jb)*z_w_con_c_full(icidx(je, jb, 1), jk, icblk(je, jb, 1)) + p_int%c_lin_e(je, 2, jb) &
    & *z_w_con_c_full(icidx(je, jb, 2), jk, icblk(je, jb, 2)))*(p_diag%vn_ie(je, jk, jb) - p_diag%vn_ie(je, jk + 1, jb)) /  &
    & p_metrics%ddqz_z_full_e(je, jk, jb))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U (no block) |
|-------|---------------------|
| `p_diag%ddt_vn_apc_pc` | `h:S  v:S  -:C` |
| `z_kin_hor_e` | `h:S  v:S` |
| `p_metrics%coeff_gradekin` | `h:S  -:C` |
| `z_ekinh` | `h:U  v:S  h:U` |
| `p_diag%vt` | `h:S  v:S` |
| `p_patch%edges%f_e` | `h:S` |
| `zeta` | `h:U  v:S  h:U` |
| `p_int%c_lin_e` | `h:S  -:C` |
| `z_w_con_c_full` | `h:U  v:S  h:U` |
| `p_diag%vn_ie` | `h:S  v:S` |
| `p_metrics%ddqz_z_full_e` | `h:S  v:S` |

</details>

### `h` full_horiz — collapsed `h:S` (2x)

6 unique arrays.

```fortran
! inside DO jb = i_startblk:i_endblk
DO je=i_startidx,i_endidx
  p_diag%vn_ie(je, 1, jb) = p_prog%vn(je, 1, jb)
  z_vt_ie(je, 1, jb) = p_diag%vt(je, 1, jb)
  z_kin_hor_e(je, 1, jb) = 0.5_wp*(p_prog%vn(je, 1, jb)**2 + p_diag%vt(je, 1, jb)**2)
  p_diag%vn_ie(je, nlevp1, jb) = p_metrics%wgtfacq_e(je, 1, jb)*p_prog%vn(je, nlev, jb) + p_metrics%wgtfacq_e(je, 2, jb) &
  & *p_prog%vn(je, nlev - 1, jb) + p_metrics%wgtfacq_e(je, 3, jb)*p_prog%vn(je, nlev - 2, jb)
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U (no block) |
|-------|---------------------|
| `p_diag%vn_ie` | `h:S  -:C` |
| `p_prog%vn` | `h:S  -:C` |
| `z_vt_ie` | `h:S  -:C` |
| `p_diag%vt` | `h:S  -:C` |
| `z_kin_hor_e` | `h:S  -:C` |
| `p_metrics%wgtfacq_e` | `h:S  -:C` |

</details>

### `v.h` full_vert, full_horiz — collapsed `h:S  v:S` (2x)

2 unique arrays.

```fortran
DO jk=1,nlev
  DO jc=i_startidx,i_endidx
    z_w_con_c(jc, jk) = p_prog%w(jc, jk, jb)
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U (no block) |
|-------|---------------------|
| `z_w_con_c` | `h:S  v:S` |
| `p_prog%w` | `h:S  v:S` |

</details>

### `h` full_horiz `accumulate` — collapsed `h:S` (1x)

7 unique arrays.

```fortran
! inside DO jb = i_startblk:i_endblk
DO je=i_startidx,i_endidx
  p_diag%vn_ie(je, 1, jb) = p_diag%vn_ie_ubc(je, 1, jb) + dt_linintp_ubc*p_diag%vn_ie_ubc(je, 2, jb)
  z_vt_ie(je, 1, jb) = p_diag%vt(je, 1, jb)
  z_kin_hor_e(je, 1, jb) = 0.5_wp*(p_prog%vn(je, 1, jb)**2 + p_diag%vt(je, 1, jb)**2)
  p_diag%vn_ie(je, nlevp1, jb) = p_metrics%wgtfacq_e(je, 1, jb)*p_prog%vn(je, nlev, jb) + p_metrics%wgtfacq_e(je, 2, jb) &
  & *p_prog%vn(je, nlev - 1, jb) + p_metrics%wgtfacq_e(je, 3, jb)*p_prog%vn(je, nlev - 2, jb)
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U (no block) |
|-------|---------------------|
| `p_diag%vn_ie` | `h:S  -:C` |
| `p_diag%vn_ie_ubc` | `h:S  -:C` |
| `z_vt_ie` | `h:S  -:C` |
| `p_diag%vt` | `h:S  -:C` |
| `z_kin_hor_e` | `h:S  -:C` |
| `p_prog%vn` | `h:S  -:C` |
| `p_metrics%wgtfacq_e` | `h:S  -:C` |

</details>

### `v.h` partial_vert, full_horiz — collapsed `h:U  v:S` (1x)

3 unique arrays.

```fortran
DO jk=nflatlev_jg,nlev
  DO jc=i_startidx,i_endidx
    z_w_concorr_mc(jc, jk) = p_int%e_bln_c_s(jc, 1, jb)*z_w_concorr_me(ieidx(jc, jb, 1), jk, ieblk(jc, jb, 1)) +  &
    & p_int%e_bln_c_s(jc, 2, jb)*z_w_concorr_me(ieidx(jc, jb, 2), jk, ieblk(jc, jb, 2)) + p_int%e_bln_c_s(jc, 3, jb) &
    & *z_w_concorr_me(ieidx(jc, jb, 3), jk, ieblk(jc, jb, 3))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U (no block) |
|-------|---------------------|
| `z_w_concorr_mc` | `h:S  v:S` |
| `p_int%e_bln_c_s` | `h:S  -:C` |
| `z_w_concorr_me` | `h:U  v:S  h:U` |

</details>

### `v.h` partial_vert, full_horiz `accumulate` — collapsed `h:S  v:S` (1x)

2 unique arrays.

```fortran
DO jk=nlev,nflatlev_jg + 1,-1
  DO jc=i_startidx,i_endidx
    z_w_con_c(jc, jk) = z_w_con_c(jc, jk) - p_diag%w_concorr_c(jc, jk, jb)
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U (no block) |
|-------|---------------------|
| `z_w_con_c` | `h:S  v:S` |
| `p_diag%w_concorr_c` | `h:S  v:S` |

</details>

### `v.h` partial_vert, full_horiz `accumulate` — collapsed `h:U  v:S` (1x)

3 unique arrays.

```fortran
DO jk=2,nlev
  DO jc=i_startidx_2,i_endidx_2
    p_diag%ddt_w_adv_pc(jc, jk, jb, ntnd) = p_diag%ddt_w_adv_pc(jc, jk, jb, ntnd) + p_int%e_bln_c_s(jc, 1, jb) &
    & *z_v_grad_w(ieidx(jc, jb, 1), jk, ieblk(jc, jb, 1)) + p_int%e_bln_c_s(jc, 2, jb)*z_v_grad_w(ieidx(jc, jb, 2), jk, ieblk(jc, &
    &  jb, 2)) + p_int%e_bln_c_s(jc, 3, jb)*z_v_grad_w(ieidx(jc, jb, 3), jk, ieblk(jc, jb, 3))
  END DO
END DO
```

<details><summary>Array access table</summary>

| Array | Role:S/U (no block) |
|-------|---------------------|
| `p_diag%ddt_w_adv_pc` | `h:S  v:S  -:C` |
| `p_int%e_bln_c_s` | `h:S  -:C` |
| `z_v_grad_w` | `h:U  v:S  h:U` |

</details>

## Widest Nest (most unique array references)

**11 unique arrays** accessed in a single nest.

- Shape: `v.h` (single-block: `v.h`)
- Collapsed: `h:U  v:S`

```fortran
DO jk=1,nlev
  DO je=i_startidx,i_endidx
    p_diag%ddt_vn_apc_pc(je, jk, jb, ntnd) = -(z_kin_hor_e(je, jk, jb)*(p_metrics%coeff_gradekin(je, 1, jb) -  &
    & p_metrics%coeff_gradekin(je, 2, jb)) + p_metrics%coeff_gradekin(je, 2, jb)*z_ekinh(icidx(je, jb, 2), jk, icblk(je, jb, 2))  &
    & - p_metrics%coeff_gradekin(je, 1, jb)*z_ekinh(icidx(je, jb, 1), jk, icblk(je, jb, 1)) + p_diag%vt(je, jk, jb) &
    & *(p_patch%edges%f_e(je, jb) + 0.5_vp*(zeta(ividx(je, jb, 1), jk, ivblk(je, jb, 1)) + zeta(ividx(je, jb, 2), jk, ivblk(je,  &
    & jb, 2)))) + (p_int%c_lin_e(je, 1, jb)*z_w_con_c_full(icidx(je, jb, 1), jk, icblk(je, jb, 1)) + p_int%c_lin_e(je, 2, jb) &
    & *z_w_con_c_full(icidx(je, jb, 2), jk, icblk(je, jb, 2)))*(p_diag%vn_ie(je, jk, jb) - p_diag%vn_ie(je, jk + 1, jb)) /  &
    & p_metrics%ddqz_z_full_e(je, jk, jb))
  END DO
END DO
```

| Array | Role:S/U (full) | Role:S/U (no block) |
|-------|-----------------|---------------------|
| `p_diag%ddt_vn_apc_pc` | `h:S  v:S  b:S  -:C` | `h:S  v:S  -:C` |
| `z_kin_hor_e` | `h:S  v:S  b:S` | `h:S  v:S` |
| `p_metrics%coeff_gradekin` | `h:S  -:C  b:S` | `h:S  -:C` |
| `z_ekinh` | `h:U  v:S  h:U` | `h:U  v:S  h:U` |
| `p_diag%vt` | `h:S  v:S  b:S` | `h:S  v:S` |
| `p_patch%edges%f_e` | `h:S  b:S` | `h:S` |
| `zeta` | `h:U  v:S  h:U` | `h:U  v:S  h:U` |
| `p_int%c_lin_e` | `h:S  -:C  b:S` | `h:S  -:C` |
| `z_w_con_c_full` | `h:U  v:S  h:U` | `h:U  v:S  h:U` |
| `p_diag%vn_ie` | `h:S  v:S  b:S` | `h:S  v:S` |
| `p_metrics%ddqz_z_full_e` | `h:S  v:S  b:S` | `h:S  v:S` |

## Array Groups (single-block)

### FULLY STRUCTURED (9 arrays)

Patterns: `h:S  v:S`

- `p_diag%w_concorr_c`
- `p_metrics%coeff1_dwdz`
- `p_metrics%coeff2_dwdz`
- `p_metrics%ddqz_z_full_e`
- `p_metrics%ddxn_z_full`
- `p_metrics%ddxt_z_full`
- `p_metrics%wgtfac_c`
- `p_metrics%wgtfac_e`
- `z_w_concorr_mc`

### FULLY STRUCTURED (5 arrays)

Patterns: `h:S  -:C`

- `p_diag%vn_ie_ubc`
- `p_int%c_lin_e`
- `p_int%e_bln_c_s`
- `p_metrics%coeff_gradekin`
- `p_metrics%wgtfacq_e`

### HAS UNSTRUCTURED (5 arrays)

Patterns: `h:S  v:S`, `h:U  v:S  h:U`

- `p_prog%w`
- `z_ekinh`
- `z_v_grad_w`
- `z_w_con_c_full`
- `z_w_concorr_me`

### FULLY STRUCTURED (4 arrays)

Patterns: `h:S  -:C`, `h:S  v:S`

- `p_diag%vn_ie`
- `p_diag%vt`
- `z_vt_ie`
- `z_w_con_c`

### FULLY STRUCTURED (4 arrays)

Patterns: `h:S`

- `p_patch%edges%f_e`
- `p_patch%edges%inv_dual_edge_length`
- `p_patch%edges%inv_primal_edge_length`
- `p_patch%edges%tangent_orientation`

### HAS UNSTRUCTURED (2 arrays)

Patterns: `h:S  -:C`, `h:S  v:S`, `h:U  v:S  h:U`

- `p_prog%vn`
- `z_kin_hor_e`

### HAS UNSTRUCTURED (2 arrays)

Patterns: `h:U  v:S  h:U`

- `z_w_v`
- `zeta`

### FULLY STRUCTURED (2 arrays)

Patterns: `h:S  v:S  -:C`

- `p_diag%ddt_vn_apc_pc`
- `p_diag%ddt_w_adv_pc`

### FULLY STRUCTURED (1 arrays)

Patterns: `-:C  h:S`

- `p_int%rbf_vec_coeff_e`

## Per-Array Summary (single-block)

| Array | Role:S/U Patterns | Conflict |
|-------|-------------------|----------|
| `p_diag%ddt_vn_apc_pc` | `h:S  v:S  -:C` |  |
| `p_diag%ddt_w_adv_pc` | `h:S  v:S  -:C` |  |
| `p_diag%vn_ie` | `h:S  -:C` \| `h:S  v:S` |  |
| `p_diag%vn_ie_ubc` | `h:S  -:C` |  |
| `p_diag%vt` | `h:S  -:C` \| `h:S  v:S` |  |
| `p_diag%w_concorr_c` | `h:S  v:S` |  |
| `p_int%c_lin_e` | `h:S  -:C` |  |
| `p_int%e_bln_c_s` | `h:S  -:C` |  |
| `p_int%rbf_vec_coeff_e` | `-:C  h:S` |  |
| `p_metrics%coeff1_dwdz` | `h:S  v:S` |  |
| `p_metrics%coeff2_dwdz` | `h:S  v:S` |  |
| `p_metrics%coeff_gradekin` | `h:S  -:C` |  |
| `p_metrics%ddqz_z_full_e` | `h:S  v:S` |  |
| `p_metrics%ddxn_z_full` | `h:S  v:S` |  |
| `p_metrics%ddxt_z_full` | `h:S  v:S` |  |
| `p_metrics%wgtfac_c` | `h:S  v:S` |  |
| `p_metrics%wgtfac_e` | `h:S  v:S` |  |
| `p_metrics%wgtfacq_e` | `h:S  -:C` |  |
| `p_patch%edges%f_e` | `h:S` |  |
| `p_patch%edges%inv_dual_edge_length` | `h:S` |  |
| `p_patch%edges%inv_primal_edge_length` | `h:S` |  |
| `p_patch%edges%tangent_orientation` | `h:S` |  |
| `p_prog%vn` | `h:S  -:C` \| `h:S  v:S` \| `h:U  v:S  h:U` |  |
| `p_prog%w` | `h:S  v:S` \| `h:U  v:S  h:U` |  |
| `z_ekinh` | `h:S  v:S` \| `h:U  v:S  h:U` |  |
| `z_kin_hor_e` | `h:S  -:C` \| `h:S  v:S` \| `h:U  v:S  h:U` |  |
| `z_v_grad_w` | `h:S  v:S` \| `h:U  v:S  h:U` |  |
| `z_vt_ie` | `h:S  -:C` \| `h:S  v:S` |  |
| `z_w_con_c` | `h:S  -:C` \| `h:S  v:S` |  |
| `z_w_con_c_full` | `h:S  v:S` \| `h:U  v:S  h:U` |  |
| `z_w_concorr_mc` | `h:S  v:S` |  |
| `z_w_concorr_me` | `h:S  v:S` \| `h:U  v:S  h:U` |  |
| `z_w_v` | `h:U  v:S  h:U` |  |
| `zeta` | `h:U  v:S  h:U` |  |
