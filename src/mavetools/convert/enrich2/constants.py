# TODO compare constants to constants in MaveCore, may be able to replace
import re

MAX_ERROR_VARIANTS = 5

supported_programs = ("enrich", "enrich2", "empiric")
extra_na = (
    "None",
    "none",
    "NONE",
    "undefined",
    "Undefined",
    "UNDEFINED",
    "na",
    "Na",
    "N/a",
    "Null",
    "",
    " ",
)
null_value_re = re.compile(r"\s+|nan|na|none|undefined|n/a|null")
surrounding_brackets_re = re.compile(r"\((.*)\)")

# HGVSP constants
hgvsp_nt_pos = "position"
hgvsp_pro_pos = "position"
hgvsp_nt_ref = "ref"
hgvsp_nt_alt = "alt"
hgvsp_pro_ref = "pre"
hgvsp_pro_alt = "post"
hgvsp_silent = "silent"


# Enrich2 constants
enrich2_synonymous = "_sy"
enrich2_wildtype = "_wt"
special_variants = ("_wt", "_sy")
synonymous_table = "synonymous"
variants_table = "variants"


# MaveDB constants
nt_variant_col = "hgvs_nt"
pro_variant_col = "hgvs_pro"
variant_columns = (nt_variant_col, pro_variant_col)
score_type = "scores"
count_type = "counts"
types = (score_type, count_type)
mavedb_score_column = "score"
