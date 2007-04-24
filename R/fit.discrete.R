`fit.discrete` <-
function(tip.data, phy, lambda=FALSE, delta=FALSE, start=1.0)
{
	if(!lambda & !delta) {
		f<-function(x) {
			likelihood.discrete(phy, tip.data, exp(x))
		}
		out<-nlm(f, p=log(start), gradtol = 1e-10)
		return(list(lnl=-out$minimum, q=exp(out$estimate), gradient=out$gradient, code=out$code, iterations=out$iterations))
	}
	
	if(lambda+delta>1) {
		cat("Cannot handle lambda and delta at the same time\n");
		return()
	}
	
	if(lambda) {
		f<-function(x) {
			likelihood.discrete(phy, tip.data, exp(x[1]), lambda=exp(x[2]))
			}
			out<-nlm(f, p=rep(log(start), 2))
			return(list(lnl=-out$minimum, q=exp(out$estimate[1]), lambda=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations))

			}
	if(delta) {
		f<-function(x) {
			likelihood.discrete(phy, tip.data, exp(x[1]), delta=exp(x[2]))
			}
			out<-nlm(f, p=rep(log(start), 2))
			return(list(lnl=-out$minimum, q=exp(out$estimate[1]), delta=exp(out$estimate[2]), gradient=out$gradient, code=out$code, iterations=out$iterations))

			}

}

