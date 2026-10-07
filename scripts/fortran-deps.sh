#!/bin/bash
#
# fortran-deps.sh - emit a make dependency rule for one Fortran object file
#
# Part of the parallel-safe build; see "Automatic Fortran module
# dependencies" in build/Makefile.
#
# Usage:
#   fortran-deps.sh STEM OUTFILE SRCDIRS "MODMAP" < source
#
#   STEM     object stem, e.g. force           -> rule "force.o: ..."
#   OUTFILE  dependency file to write, e.g. force.d (updated atomically)
#   SRCDIRS  source directories (the Makefile VPATH); used to recognise
#            project modules. A used module with no source on disk is
#            assumed external (mpi, omp_lib, hdf5, ...) and gets no edge.
#   MODMAP   space-separated "mod=obj[,obj...]" pairs. Wins over the fixed
#            table and the same-name default below. Carries the per-SETUP
#            families where only one variant is compiled at a time
#            (setup, bondiexact, analysis, moddump, kernel, metric,
#            inject, externalforces, extern_binary, testexternf).
#
# Reads preprocessed (or raw) Fortran source on stdin. Extracts "use MOD"
# edges and "submodule (PARENT)" edges; drops "use, intrinsic",
# self-edges and externals. Unknown modules with no source on disk print
# a stderr warning (add them to MODMAP or the FIX table if project code).
#
# Mapping order (first hit wins):
#   1. dynamic MODMAP from the Makefile (setup-dependent families)
#   2. FIX table below (renames; generated from the source tree, refresh
#      if modules are renamed, e.g. densityforce=dens, dim=config)
#   3. same-name default (module foo -> foo.o), kept only if a matching
#      source exists on disk, else assumed external (quietly for known
#      compiler/library modules, warning otherwise).
#
set -u

STEM="$1"
OUT="$2"
SRCDIRS="$3"
MODMAP="$4"

TMP="$OUT.$$"
trap 'rm -f "$TMP"' EXIT

# lowercase stems of all Fortran sources on disk (single pass, no per-file fork)
DISK=""
for _d in $SRCDIRS; do
    _dd=$(echo "$_d" | tr -d '"')
    [ -n "$_dd" ] && [ -d "$_dd" ] || continue
    DISK="$DISK $(ls -1 "$_dd" 2>/dev/null | sed -n 's/\.[fF]90$//p' | tr 'A-Z' 'a-z')"
done

awk -v stem="$STEM" -v modmap="$MODMAP" -v disk="$DISK" '
function warn(msg) { print "fortran-deps.sh: warning: " msg > "/dev/stderr" }
function adddep(mod,   objs, n, arr, j, obj) {
    if (mod == SELF) return
    if (mod in DYN) objs = DYN[mod]
    else if (mod in FIX) objs = FIX[mod]
    else if (mod in QUIETEXT) return
    else if (mod in ONDISK) objs = mod
    else { warn(stem ": no source for module \"" mod "\", assuming external"); return }
    n = split(objs, arr, ",")
    for (j = 1; j <= n; j++) {
        obj = arr[j]
        if (obj == "") continue
        if (tolower(obj) == SELF) continue
        if (!(obj in SEEN)) { SEEN[obj] = 1; DEPS[++NDEP] = obj }
    }
}
function doprocess(line,   l, s) {
    # OpenMP conditional-compilation sentinel ("!$ use dim") is code
    # under -fopenmp; keep it (an extra edge is harmless without OpenMP)
    sub(/^[ \t]*!\$[ \t]/, "", line)
    sub(/!.*$/, "", line)
    l = tolower(line)
    if (match(l, /^[ \t]*submodule[ \t]*\([ \t]*[a-z_][a-z_0-9]*/)) {
        s = substr(l, RSTART, RLENGTH)
        sub(/.*\(/, "", s); gsub(/[ \t]/, "", s)
        adddep(s)
    }
    if (l ~ /^[ \t]*use([ \t,]|$)/) {
        sub(/^[ \t]*use[ \t]*/, "", l)
        sub(/^,[ \t]*non_intrinsic[ \t]*::[ \t]*/, "", l)
        if (l ~ /^,[ \t]*intrinsic[ \t]*::/) return
        if (match(l, /^[a-z_][a-z_0-9]*/)) adddep(substr(l, RSTART, RLENGTH))
    }
}
BEGIN {
    SELF = tolower(stem)
    FIX["allocutils"]="utils_allocate"
    FIX["amusephantom"]="libphantom-amuse"
    FIX["boundary_dyn"]="boundary_dynamic"
    FIX["chem"]="h2chem"
    FIX["cooling_gammie_pl"]="cooling_gammie_PL"
    FIX["cpuinfo"]="utils_cpuinfo"
    FIX["cubic"]="cubicsolve"
    FIX["datautils"]="utils_datafiles"
    FIX["densityforce"]="dens"
    FIX["derivutils"]="utils_deriv"
    FIX["dim"]="config"
    FIX["discanalysisutils"]="utils_disc"
    FIX["dtypekdtree"]="dtype_kdtree"
    FIX["dump_utils"]="utils_dumpfiles"
    FIX["easter_egg"]="egg"
    FIX["eos_hiir"]="eos_HIIR"
    FIX["ephemeris"]="utils_ephemeris"
    FIX["ev2dotutils"]="ev2dot_utils"
    FIX["evolveplanet"]="evolve_planet"
    FIX["evutils"]="utils_evfiles"
    FIX["extern_bfield"]="extern_Bfield"
    FIX["fileutils"]="utils_filenames"
    FIX["forces"]="force"
    FIX["getneighbours"]="utils_getneighbours"
    FIX["gravwaveutils"]="utils_gravwave"
    FIX["healpix"]="utils_healpix"
    FIX["hiiregion"]="H2regions"
    FIX["implicit"]="utils_implicit"
    FIX["infile_utils"]="utils_infiles"
    FIX["injectutils"]="utils_inject"
    FIX["interpolations3d"]="interpolate3D"
    FIX["interpolations3d_amr"]="interpolate3D_amr"
    FIX["io_summary"]="utils_summary"
    FIX["ionization_mod"]="ionization"
    FIX["krome_interface"]="krome"
    FIX["libphantomsplash"]="libphantom-splash"
    FIX["linalg"]="utils_linalg"
    FIX["mathfunc"]="utils_mathfunc"
    FIX["mesa_microphysics"]="eos_mesa_microphysics"
    FIX["metric_interp"]="interp_metric"
    FIX["mpc"]="utils_mpc"
    FIX["mpibalance"]="mpi_balance"
    FIX["mpidens"]="mpi_dens"
    FIX["mpiderivs"]="mpi_derivs"
    FIX["mpidomain"]="mpi_domain"
    FIX["mpiforce"]="mpi_force"
    FIX["mpimemory"]="mpi_memory"
    FIX["mpitree"]="mpi_tree"
    FIX["mpiutils"]="mpi_utils"
    FIX["neighkdtree"]="neigh_kdtree"
    FIX["nicil_sup"]="nicil_supplement"
    FIX["omputils"]="utils_omp"
    FIX["orbits"]="utils_orbits"
    FIX["raytracer"]="utils_raytracer"
    FIX["raytracer_all"]="utils_raytracer_all"
    FIX["relaxstar"]="relax_star"
    FIX["rho_profile"]="density_profiles"
    FIX["setbfield"]="set_Bfield"
    FIX["setbinary"]="set_binary"
    FIX["setcubiccore"]="set_cubic_core"
    FIX["setdisc"]="set_disc"
    FIX["setfixedentropycore"]="set_fixedentropycore"
    FIX["sethier_utils"]="set_hierarchical_utils"
    FIX["sethierarchical"]="set_hierarchical"
    FIX["setorbit"]="set_orbit"
    FIX["setplanets"]="set_planets"
    FIX["setplummer"]="set_plummer"
    FIX["setshock"]="set_shock"
    FIX["setsoftenedcore"]="set_softened_core"
    FIX["setsolarsystem"]="set_solarsystem"
    FIX["setstar"]="set_star"
    FIX["setstar_utils"]="set_star_utils"
    FIX["setunits"]="set_units"
    FIX["setvfield"]="set_vfield"
    FIX["slab"]="set_slab"
    FIX["sortutils"]="utils_sort"
    FIX["spherical"]="set_sphere"
    FIX["sphngutils"]="utils_sphNG"
    FIX["splineutils"]="utils_spline"
    FIX["splitmergeutils"]="utils_splitmerge"
    FIX["step_lf_global"]="step_leapfrog"
    FIX["structurefn_part"]="struct_part"
    FIX["systemutils"]="utils_system"
    FIX["table_utils"]="utils_tables"
    FIX["test"]="testsuite"
    FIX["testapr"]="test_apr"
    FIX["testbinary"]="test_binary"
    FIX["testcoala"]="test_coala"
    FIX["testcooling"]="test_cooling"
    FIX["testcorotate"]="test_corotate"
    FIX["testdamping"]="test_damping"
    FIX["testderivs"]="test_derivs"
    FIX["testdust"]="test_dust"
    FIX["testeos"]="test_eos"
    FIX["testeos_stratified"]="test_eos_stratified"
    FIX["testgeometry"]="test_geometry"
    FIX["testgnewton"]="test_gnewton"
    FIX["testgr"]="test_gr"
    FIX["testgravity"]="test_gravity"
    FIX["testgrowth"]="test_growth"
    FIX["testindtstep"]="test_indtstep"
    FIX["testiorig"]="test_iorig"
    FIX["testkdtree"]="test_kdtree"
    FIX["testkernel"]="test_kernel"
    FIX["testlinalg"]="test_linalg"
    FIX["testlum"]="test_luminosity"
    FIX["testmpi"]="test_mpi"
    FIX["testneigh"]="test_neigh"
    FIX["testnimhd"]="test_nonidealmhd"
    FIX["testorbits"]="test_orbits"
    FIX["testpart"]="test_part"
    FIX["testpoly"]="test_poly"
    FIX["testptmass"]="test_ptmass"
    FIX["testradiation"]="test_radiation"
    FIX["testrwdump"]="test_rwdump"
    FIX["testsedov"]="test_sedov"
    FIX["testsetdisc"]="test_setdisc"
    FIX["testsethier"]="test_hierarchical"
    FIX["testsetstar"]="test_setstar"
    FIX["teststep"]="test_step"
    FIX["testunits"]="test_units"
    FIX["testutils"]="utils_testsuite"
    FIX["testwind"]="test_wind"
    FIX["timestep_ind"]="utils_indtimesteps"
    FIX["timing"]="utils_timing"
    FIX["unifdis"]="set_unifdis"
    FIX["vectorutils"]="utils_vectors"
    FIX["velfield"]="velfield_fromcubes"
    QUIETEXT["mpi"]=1; QUIETEXT["mpi_f08"]=1; QUIETEXT["omp_lib"]=1
    QUIETEXT["iso_c_binding"]=1; QUIETEXT["iso_fortran_env"]=1
    QUIETEXT["ieee_arithmetic"]=1; QUIETEXT["ieee_exceptions"]=1
    QUIETEXT["hdf5"]=1; QUIETEXT["giza"]=1; QUIETEXT["mcfost2phantom"]=1
    QUIETEXT["krome_user"]=1; QUIETEXT["krome_main"]=1
    n = split(modmap, m, " ")
    for (i = 1; i <= n; i++) {
        if (m[i] == "") continue
        e = index(m[i], "=")
        if (e == 0) continue
        DYN[tolower(substr(m[i], 1, e-1))] = substr(m[i], e+1)
    }
    for (k in DYN) if (DYN[k] == "") delete DYN[k]
    nd = split(disk, dd, " ")
    for (i = 1; i <= nd; i++) if (dd[i] != "") ONDISK[dd[i]] = 1
    BUF = ""
    NDEP = 0
}
{
    line = $0
    if (BUF != "") { sub(/^[ \t]*&/, "", line); line = BUF line; BUF = "" }
    if (line ~ /&[ \t]*$/) { sub(/&[ \t]*$/, "", line); BUF = line; next }
    doprocess(line)
}
END {
    if (BUF != "") doprocess(BUF)
    printf "%s.o:", stem
    for (i = 1; i <= NDEP; i++) printf " %s.o", DEPS[i]
    printf "\n"
}
' > "$TMP" && mv "$TMP" "$OUT"
