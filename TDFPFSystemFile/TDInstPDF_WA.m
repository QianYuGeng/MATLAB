function [Xfail_dv,fxhx_no,Eu,Betaint,Ifail] = TDInstPDF_WA(Probs,N,Method)
Ifail_old=false(N,1);
Xfail_dv=cell(Probs.Nt,1);
fxhx_no=cell(Probs.Nt,1);
for ind_t=1:Probs.Nt
[Xfail_dv{ind_t,1},fxhx_no{ind_t,1},Eu,Betaint,Ifail_old] = InstrumentalPDF_WA_t(Probs,N,Method,ind_t,Ifail_old);
Ifail(:,ind_t)=Ifail_old;
end 
end