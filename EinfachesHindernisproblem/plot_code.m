figure(2);
subplot(2,1,1);
plot(1:recursion_depth+3,osc1_term,':o',1:recursion_depth+3,osc2_term, '-.*',1:recursion_depth+3,osc_term,':x');
ymin = min([min(osc_term_uniform),min(osc1_term_uniform),min(osc2_term_uniform)])-5;
ymax = max([max(osc_term_uniform),max(osc1_term_uniform),max(osc2_term_uniform)])+5;
axis([0.5,recursion_depth+3.5,ymin,ymax]);
legend('osc1','osc2','oscillation','location','best');

subplot(2,1,2);
plot(1:recursion_depth+3,J_error,'--o',...
    1:recursion_depth+3,IQ_plot,'-.*',1:recursion_depth+3,rhoS_plot,'-.x');
ymin = min([min(J_error),min(IQ_plot),min(rhoS_plot)])-0.5;
ymax = max([max(J_error),max(IQ_plot),max(rhoS_plot)])+0.5;
axis([0.5,recursion_depth+3.5,ymin,ymax]);
legend('functional error','estimated error','error indicator','location',...
    'best');