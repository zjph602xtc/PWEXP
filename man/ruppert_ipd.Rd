\name{ruppert_ipd}
\alias{ruppert_ipd}
\title{OS Curves from Ruppert (2020)}
\usage{ruppert_ipd
       data(ruppert_ipd, package="PwePred")
}
\description{OS Curves for IPI 4-5 (red curve) from Figure 2A in Ruppert (2020)}
\format{
\tabular{ll}{
Survival.time:\tab time\cr
Status:\tab 1 for events, 0 for censoring
}}
\note{We digitize the OS curves for IPI 4-5 (red curve) from Figure 2A in Ruppert (2020),  using WebPlotDigitizer (Rohatgi A (2024)). Subsequently, we reconstruct individual subject data (IPD) using the \pkg{IPDfromKM} package (Liu (2021)).
   }
\references{
Ruppert AS, Dixon JG, Salles G, Wall A, Cunningham D, Poeschel V, Haioun C, Tilly H, Ghesquieres H, Ziepert M, et al. (2020). "International Prognostic Indices in Diffuse Large B-cell Lymphoma: A Comparison of IPI, R-IPI, and NCCN-IPI." Blood, The Journal of
the American Society of Hematology, 135(23), 2041-2048.

Rohatgi A (2024). "Webplotdigitizer: Version 4.7." URL: https://automeris.io/WebPlotDigitizer.

Liu N, Zhou Y, Lee JJ (2021). "IPDfromKM: Reconstruct Individual Patient Data from Published Kaplan-Meier Survival Curves." BMC medical Research Methodology, 21(1), 111.
  }
\keyword{PwePred}
