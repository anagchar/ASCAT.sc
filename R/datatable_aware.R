## The package uses data.table syntax (:=, .I, .N) but does not import data.table in NAMESPACE.
## Without this flag, data.table treats calls from the ASCAT.sc namespace as not data.table-aware
## (cedta() is FALSE) and `dt[, x := ...]` fails, e.g. in ascatsc_plot().
.datatable.aware <- TRUE
