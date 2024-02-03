# @name Region
# @export
Region = setRefClass("Region")

# @name UnivariateConstRegion
# @export
UnivariateConstRegion = setRefClass("UnivariateConstRegion",
	contains = "Region",
	fields = c("a", "b", "w", "g", "log_w_max", "log_w_min", "log_prob"))

# @name IntUnivariateConstRegion
# @export
IntUnivariateConstRegion = setRefClass("IntUnivariateConstRegion",
	contains = "Region",
	fields = c("a", "b", "w", "g", "log_w_max", "log_w_min", "log_prob"))

# @exportClass FMMProposal
FMMProposal = setRefClass("FMMProposal",
	fields = c("regions", "log_xi_upper", "log_xi_lower", "bifurcatable"))
