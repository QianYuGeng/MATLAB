function [Pft, CovPft] = TDPFsita_WA_BOC(Sita,Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,flag)
global Prob indt
for indt=1:Prob.Nt
[Pft(:,indt), CovPft(:,indt)] = PFsita_WA_BOC(Sita,Xfail(indt,:),fxhx_no(indt,:),Eu,Betaint,Probk,N,Method,flag);
end
end