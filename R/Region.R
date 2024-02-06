#' Region
#' @description
#' A Region contains all of the problem-specific logic for the vws sampler.
#' This is an "abstract" R6 class that defines the interface for a Region.
Region = R6::R6Class(classname = "Region", portable = TRUE, lock_class = TRUE,
public = list(

#' @description
#' Density function \eqn{g} for the base distribution.
#'
#' @param x Density argument.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
d_base = function(x, log = FALSE)
{
	stop("d_base: abstract interface method")
},

#' @description
#' Weight function \eqn{w}.
#'
#' @param x Argument to weight function.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
w = function(x, log = TRUE)
{
	stop("w: abstract interface method")
},

#' @description
#' Generate a draw from \eqn{g_j} specific to this region.
#'
#' @param n Number of draws to generate.
#'
#' @return A list of draws, with one draw per list element.
r = function(n)
{
	stop("r: abstract interface method")
},

#' @description
#' Density of \eqn{g_j} specific to this region.
#'
#' @param x Density argument.
d = function(x)
{
	stop("d: abstract interface method")
},

#' @description
#' Test if given \code{x} is in the support for the \eqn{g_j} specific to this
#' region.
#'
#' @param x Density argument.
s = function(x)
{
	stop("s: abstract interface method")
},

#' @description
#' Majorized weight function \eqn{\overline{w}_j} for this region.
#'
#' @param x Argument to weight function.
#' @param log logical; if \code{TRUE}, return result on the log-scale.
w_major = function(x, log = TRUE)
{
	stop("w_major: abstract interface method")
},

#' @description
#' Bifurcate this region into two regions. Use \code{x} as the bifurcation
#' point if it is not \code{NULL}. Otherwise, select a point for bifurcation.
#'
#' @param x An optional bifurcation point.
bifurcate = function(x = NULL)
{
	stop("bifurcate: abstract interface method")
},

#' @description
#' The quantity \eqn{\overline{\xi}_j} for this region.
#'
#' @param log logical; if \code{TRUE}, return result on the log-scale.
xi_upper = function(log = TRUE)
{
	stop("xi_upper: abstract interface method")
},

#' @description
#' The quantity \eqn{\underline{\xi}_j} for this region.
#'
#' @param log logical; if \code{TRUE}, return result on the log-scale.
xi_lower = function(log = TRUE)
{
	stop("xi_lower: abstract interface method")
}

) # Close public
) # Close class

