"""A clean environment for the shell / perl / R based tools (IEDB bundles, MixMHCpred, PRIME).

A shell with conda active puts conda's libraries first (LD_LIBRARY_PATH) and its perl / R settings in the environment. The system programs the
tools call (tcsh, gawk, perl, R) then load conda's libraries instead of their own: tcsh prints "libtinfo.so.6: no version information available",
and perl or R can fail silently, which showed up as CD4episcore returning nothing. These tools get only what they need.
"""
import os

DROP = {"PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV", "LD_LIBRARY_PATH", "LD_PRELOAD", "DYLD_LIBRARY_PATH", "PERL5LIB", "PERLLIB", "PERL_LOCAL_LIB_ROOT",
        "PERL_MB_OPT", "PERL_MM_OPT", "R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_PROFILE", "R_PROFILE_USER", "TCLLIBPATH", "PKG_CONFIG_PATH"}


def clean(path_first=None):
    e = {k: v for k, v in os.environ.items() if k not in DROP and not k.startswith(("CONDA_", "_CE_"))}
    e["PATH"] = (path_first + ":" if path_first else "") + "/usr/local/bin:/usr/bin:/bin"
    return e
