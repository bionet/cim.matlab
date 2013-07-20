classdef cim_utility
    methods (Static)
        function runtime(msg,time) 
            display([msg num2str(floor(time/60)) ''' ' ...
                      num2str(time - 60*floor(time/60), '%3.1f') '"']);
        end
        function u = synthesizeSignal(t_fft,t_sig,c,L,T,varargin)
            % t_fft is the time course of which the fft coefficient is computed
            % t_sig is the time course of the output signal
            if nargin == 6 && ~strcmpi(varargin,'normalize');
                error('Unknown optional parameter. Use "normalize" only...');
            elseif nargin >= 7
                error('Too many input arguments....');
            end
            dt = t_fft(2)-t_fft(1);
            u = zeros( size(c,1),numel(t_fft) );
            for i = 1:size(c,1)
                c_pad   = [c(i,L+1:end) zeros(1,numel(t_fft)-2*L-1) c(i,1:L)];
                u(i,:)  = real(ifft(c_pad))*numel(t_fft)/sqrt(T);
            end
            idx = floor(mod(t_sig(1)-t_fft(1),T)/dt)+1;
            u = u(:, [idx:end 1:idx-1]);
            if nargin == 6
                u = diag(1./max(abs(u),[],2))*u;     
            end
        end
        function plotFilter(t,h,Ph)
            plot(t,h,'--k',t,Ph,'r');
            xlabel('Time, [s]');ylabel('Amplitude');xlim([t(1)  t(end)]);
            title('Filter and its projection');
            legend('$h$','$\mathcal{P}h$');
            if max(h) > 1e4, 
                hold on;plot(0,1.08*max(Ph),'^k','MarkerFaceColor','k');
                ylim([-50 1.1*max(Ph)]); 
            else
                ylim([1.1*min(h) 1.1*max(h)]);
            end
        end
        function plotInputSignal( t, u)
            cs = colormap(lines(64));
            N = size(u,1);
            if N<=6
                for i=1:N-1
                    subplot(N,1,i); plot(t,u(i,:),'color',cs(i,:));
                        ylabel('Amplitude');xlim([t(1) t(end)]);
                        set(gca,'xticklabel',[]);
                        legend(sprintf('Input signal $u^%d(t)\\quad$',i));
                end
                subplot(N,1,N); plot(t,u(N,:),'color',cs(N,:));
                    ylabel('Amplitude');xlim([t(1) t(end)]);xlabel('Time, [s]');
                    legend(sprintf('Input signal $u^%d(t)\\quad$',N));
            elseif N>6 && N<=12
                for i=1:N
                    subplot(ceil(N/2),2,i); plot(t,u(i,:));
                end
            else
                for i=1:N
                    subplot(ceil(N/3),3,i); plot(t,u(i,:));
                end
            end
        end
        function plotFilterOutput( t, u, v )
            N = size(u,1);
            cs = colormap(lines(64));
            if N<=6
                for i=1:N-1
                    subplot(N,2,2*i-1); plot(t,u(i,:),'color',cs(i,:));
                        ylabel('Amplitude');xlim([t(1) t(end)]);
                        set(gca,'xticklabel',[]);
                        legend(sprintf('Input signal $u^%d(t)\\quad$',i));
                    subplot(N,2,2*i);   plot(t,v(i,:),'color',cs(i,:));
                        xlim([t(1) t(end)]);set(gca,'xticklabel',[]);
                        legend(sprintf('Filter output $v^%d(t)\\quad$',i));
                end
                subplot(N,2,2*N-1); plot(t,u(N,:),'color',cs(N,:));
                    ylabel('Amplitude');xlim([t(1) t(end)]);xlabel('Time, [s]');
                    legend(sprintf('Input signal $u^%d(t)\\quad$',N));
                subplot(N,2,2*N); plot(t,v(N,:),'color',cs(N,:));
                    xlim([t(1) t(end)]);xlabel('Time, [s]');
                    legend(sprintf('Filter output $v^%d(t)\\quad$',N));
            elseif N>6 && N<=12
                for i=1:N
                    subplot(ceil(N/2),2,i); plot(t,u(i,:));
                end
            else
                for i=1:N
                    subplot(ceil(N/3),3,i); plot(t,u(i,:));
                end
            end                
        end
        function plotEncodingResult(f,L,t,u,voltage,s,delta_v,bias,kappa,delta,sigma_delta)
            intvl  = [t(1) t(end)];
            tk_idx = s(1,:) == 1;
            subplot(311);
                plot(t,u(1,:));
                ylabel('Amplitude');xlim(intvl);set(gca,'xticklabel',[]);
                title('(a)$\qquad$Input signal $u^1(t)$');
                legend(['$\Omega = 2\pi\cdot$' num2str(f) 'rad/s, $L = ' num2str(L) '\quad$'],...
                        'location','East');
            subplot(312);
                plot(t , voltage(1,:),'-b'); 
                hold on; plot(intvl, delta*ones(1,2),'--r');
                hold on; plot(t(tk_idx), delta_v{1}(1:end-1),'or',...
                              'linewidth', 1, 'markersize',6, 'markerfacecolor','r'); 
                ylabel('Membrane Potential');set(gca,'xticklabel',[]);
                title('(b)$\qquad$Integrator output vs. Time');
                if sigma_delta > 1e-20
                    legend('$\int_{t_k}^{t}(u\ast h)(s)ds,\,\forall k\qquad$',...
                            ['$\delta = ' num2str(delta) '$'],...
                            ['$\Delta_k\sim\mathcal{N}\big(\delta, (', num2str(sigma_delta), '\delta)^2\big)\quad$'],...
                            'Location','southeast');
                else
                    legend('$\int_{t_k}^{t}(u\ast h)(s)ds,\,\forall k\qquad$',...
                           ['$\delta = ' num2str(delta) '$'],...
                           '$\int_{t_k}^{t}(u\ast h)(s)ds=\delta$',...
                           'Location','southeast');
                end
                xlim(intvl);ylim([-0.1*delta 1.2*delta]);
                
            subplot(313);
                stem(t(tk_idx), 1.1*ones(1,sum(tk_idx)), '^k', 'filled');
                box on; set(gca,'yticklabel',[],'ytick',[]);
                title('(c)$\qquad$IAF spike train for $u^1(t)$');  
                xlabel('Time, [s]'); xlim(intvl);ylim([0 1.2]);
                legend(['Spikes, $n = ' num2str(sum(s(1,:)),'%3.0f') '\quad$'],...
                       'Location','East'); 
        end
        function [Fx,freq]=centeredFFT(x,Fs,N)
            %this is a custom function that helps in plotting the two-sided spectrum
            %x is the signal that is to be transformed
            %Fs is the sampling rate
            %N is the number of points returned in the FFT result

            %this part of the code generates that frequency axis
            freq = (-floor(N/2):N-floor(N/2)-1)*Fs/N;
            Fx = fftshift(fft(x,N)/numel(x));
        end
        function plotPaperFig(T1,T2,T,f,W,L,bias,delta,sigma_delta,dt,...
                         t_filt,h,Ph,t,u,v,...
                         s,voltage,delta_v,h_rec,total_spikes_used)
            % compute the recovery error
            idx = find( t_filt >= T1 & t_filt <= T2 );
            h_hhat_err   = (h-h_rec)/max(abs(h));
            ph_hhat_err  = (Ph-h_rec)/max(abs(Ph));
            h_hhat_rmse  = sqrt( dt*trapz(h_hhat_err(idx).^2)/dt/numel(idx));
            ph_hhat_rmse = sqrt( dt*trapz(ph_hhat_err(idx).^2)/dt/numel(idx));
            % err_rms_norm_2 = sqrt(dt*trapz(err_norm_2(t_frame_i).^2)/(t_frame(end)-t_frame(1)));
            % shift time course to start at zero
            t = t - t(1);
            intvl = [t(1) t(end)];
            cs = colormap(lines(64));
            subplot(421);
                plot(t, u); 
                xlim(intvl);ylabel('Amplitude');set(gca,'xticklabel',[]);
            subplot(423);
                plot(t, v(1,:)+bias,intvl, bias*ones(1,2),'--r',intvl, zeros(1,2),'--k');
                ylabel('Amplitude');xlim(intvl);set(gca,'xticklabel',[]);
            subplot(425);
                plot(t, voltage(1,:),'-b',intvl, delta*ones(1,2),'--r');hold on;
                plot(t(logical(s(1,:))), delta_v{1}(1:end-1), 'or',...
                    'linewidth', 1, 'markersize',6, 'markerfacecolor','r');
                ylabel('Amplitude');xlim(intvl);ylim([0 1.2*delta]); set(gca,'xticklabel',[]);
            subplot(427);
                if size(u,1) == 1, cs(1,:) = [0 0 0]; end
                for i=1:size(u,1)
                    tk_idx = s(i,:)==1;
                    hold on;stem(t(tk_idx) , 1.1*ones(1,sum(tk_idx)),...
                                 '^', 'fill', 'color', cs(i,:) ); 
                end
                set(gca,'yticklabel',[],'ytick',[],'box','on');
                xlabel('Time, [s]');xlim(intvl);ylim([0 1.2]);
            subplot(422);
                if max(h) > 1e4,
                    plot([t_filt(1) t_filt(end)], zeros(1,2), '--k', ...
                         t_filt, Ph, '-b', t_filt, h_rec, '-r',...
                         [0 0],[0 max(Ph)],'--k',0,max(Ph),'^k','MarkerFaceColor','k');
                    axis([t_filt(1) t_filt(end) -50 1.1*max(Ph)]); 
                else
                    plot(t_filt, h, '--k', t_filt, Ph, '-b', t_filt, h_rec, '-r');
                    axis([t_filt(1) t_filt(end) 1.1*min(h) 1.1*max(h)]);
                end
                legend(['$h,\,$MSE$(\mathcal{P}h^*,h)$ =' num2str(10*log10(h_hhat_rmse^2), '%3.1f') 'dB $\quad$'],...
                       ['$\mathcal{P}h,\,$MSE$(\mathcal{P}h^*,\mathcal{P}h)$ =' num2str(10*log10(ph_hhat_rmse^2), '%3.1f') 'dB $\qquad$'],...
                       ['$\mathcal{P}h^*$, from $n = ' num2str(total_spikes_used) '$ spikes'],...
                       'Location','NorthEast');
                title('$(e)\qquad$Original filter vs. the identified filter');
                ylabel('Amplitude');xlabel('Time, [s]');
                
            subplot(424);
                t_fft = -11*T/2:dt:9*T/2-dt;
                K = sinc((2*L+1)*W/(2*L)*t_fft/pi)./sinc(W/(2*L)*t_fft/pi)*(2*L+1)/T;
                n_fft = 2^(1+nextpow2(numel(K)));
                [FK freq_K] = cim_utility.centeredFFT(K, 1/dt, n_fft);
                plot(freq_K, 10*log10(abs(FK)));
                title('$(f)\qquad$Fourier amplitude spectrum of $K$');grid on;
                legend('supp$(\mathcal{F}K)=[-\Omega,\,\Omega]\quad$','location','southeast');
                ylabel('$10\log|\mathcal{F}K|\;$'); set(gca,'xticklabel',[]);
                xlim([-125 125]);ylim([-40 20]);
                
            subplot(426);
                h_long = [ zeros(1,numel(find(t_fft<t_filt(1)))) h zeros(1,numel(find(t_fft>t_filt(end))))];
                [Fh freq_h] = cim_utility.centeredFFT(h_long, 1/dt, n_fft);
                if max(h) > 1e4
                    plot(freq_h, zeros(size(freq_h)));
                    legend('supp$(\mathcal{F}h)\supset[-\Omega,\,\Omega]\quad$','location','southeast');
                    
                else
                    plot(freq_h, 10*log10(abs(Fh)));
                    legend('supp$(\mathcal{F}h)\supset[-\Omega,\,\Omega]\quad$','location','southeast');
                end
                set(gca,'xticklabel',[]);grid on;
                xlim([-125 125]);ylim([-40 20]);ylabel('$10\log|\mathcal{F}h|\;$');
                title('$(g)\qquad$Fourier amplitude spectrum of $h$');

            subplot(428);
                h_rec_rep   = repmat(h_rec(1:end-1),1,10);
                [Fh_rec freq_h_rec] = cim_utility.centeredFFT(h_rec_rep, 1/dt, n_fft);
                plot(freq_h_rec, 10*log10(abs(Fh_rec)));
                xlim([-125 125]);ylim([-40 20]);grid on;
                ylabel('$10\log|\mathcal{F}\hat{h}|\;$');xlabel('Frequency, [Hz]');
                title('$(h)\qquad$Fourier amplitude spectrum of $\mathcal{P}h^*$');
                legend('supp$(\mathcal{F}\mathcal{P}h^*)=[-\Omega,\,\Omega]\quad$','location','southeast');

            if size(u,1) == 1
                subplot(421);
                title('$(a)\qquad$Input signal $u(t)$');
                legend(['$\Omega = 2\pi\cdot$' num2str(f) 'rad/s, $L = ' num2str(L) '\quad$'],...
                           'location','northeast');
                subplot(423);
                title('$(b)\qquad$Biased filter output $v(t)+b$');
                legend('$(u\ast h)(t)+b\quad$', ['Bias $b = ' num2str(bias, '%3.2f') '$'],...
                       'Zero line','location','southeast');

                subplot(425);
                title('(c)$\qquad$Ideal IAF neuron response');
                if sigma_delta > 1e-20
                     legend('$\int_{t_k}^{t}(u\ast h)(s)ds,\,\forall k\qquad$',...
                            ['$\delta = ' num2str(delta) '$'],...
                            ['$\Delta_k\sim\mathcal{N}\big(\delta, (', num2str(sigma_delta), '\delta)^2\big)\quad$'],...
                            'Location','southeast');
                else
                    legend('$\int_{t_k}^{t}(u\ast h)(s)ds,\,\forall k\qquad$',...
                           ['$\delta = ' num2str(delta) '$'],...
                           '$\int_{t_k}^{t}(u\ast h)(s)ds=\delta$',...
                           'Location','southeast');
                end

                subplot(427);
                title('$(d)\qquad$Output sequence $(t_k)_{k=1}^n$');
                legend(['Spikes, $n = ' num2str(total_spikes_used,'%3.0f') '\quad$'],...
                       'Location','East');            
            else
                subplot(421);
                    title(['$(a)\qquad$Input signals $\{u^i\}_{i=1}^{' num2str(size(u,1)) '}$. $\Omega = 2\pi\cdot$' num2str(f) 'rad/s, $L = ' num2str(L) '\quad$']);
                    if size(u,1)<=5
                        legend_text = cell(size(u,1),1);
                        for i=1:size(u,1)
                            legend_text{i} = ['$u^' num2str(i) '(t)$'];
                        end
                        legend(legend_text,'location','east');
                    else
                        legend(['$u^1(t),\;\dots\;,u^{' num2str(size(u,1)) '}(t)$' ],'location','east');
                    end
                subplot(423);
                    title('$(b)\qquad$Biased filter output $v^1(t)+b$');
                    legend('$(u^1\ast h)(t)+b\quad$', ['Bias $b = ' num2str(bias, '%3.2f') '$'], 'Zero line','location','southeast');
                subplot(425);
                title('(c)$\qquad$Ideal IAF neuron response to $u^1$');
                    if sigma_delta > 1e-20
                         legend('$\int_{t_k}^{t}(u^1\ast h)(s)ds,\,\forall k\qquad$',...
                                ['$\delta = ' num2str(delta) '$'],...
                                ['$\delta_k\sim\mathcal{N}\big(\delta, (', num2str(sigma_delta), '\delta)^2\big)\quad$'], 'Location','east');
                    else
                        legend('$\int_{t_k}^{t}(u^1\ast h)(s)ds,\,\forall k\qquad$',...
                               ['$\delta = ' num2str(delta) '$'], '$\int_{t_k}^{t}(u^1\ast h)(s)ds=\delta$', 'Location','southeast');
                    end

                subplot(427);
                title(['$(d)\qquad$Output sequences $(t^i_k)_{k=1}^{n^i}$. $\sum_{i=1}^{' num2str(size(u,1)) '}n^i = ' num2str(total_spikes_used,'%3.0f') '\quad$']);
                if size(u,1)<=5
                    legend_text = cell(size(u,1),1);
                    for i=1:size(u,1)
                        tk_num = sum(s(i,:)==1);
                        legend_text{i} = ['$(t_k^' num2str(i) ')_{k=1}^{' num2str(tk_num) '}$'];
                    end
                    legend(legend_text,'Location','East');
                else
                    legend(['$(t_k^1)_{k=1}^{' num2str(length(tk{1})) '}\;,\dots\;,\;(t_k^{' num2str(size(u_simul,1)) '})_{k=1}^{' num2str(length(tk{size(u_simul,1)})) '}\;$' ],'location','east');
                end
            end
        end
        function plotPaperFig8(T1,T2,T,f,W,L,mu,bias,dt,lc,...
                         t_filt,h,Ph,t,u,v,y1,y2,s,h_rec)
            % compute the recovery error
            idx = find( t_filt >= T1 & t_filt <= T2 );
            h_hhat_err   = (h-h_rec)/max(abs(h));
            ph_hhat_err  = (Ph-h_rec)/max(abs(Ph));
            h_hhat_rmse  = sqrt( dt*trapz(h_hhat_err(idx).^2)/dt/numel(idx));
            ph_hhat_rmse = sqrt( dt*trapz(ph_hhat_err(idx).^2)/dt/numel(idx));
            % err_rms_norm_2 = sqrt(dt*trapz(err_norm_2(t_frame_i).^2)/(t_frame(end)-t_frame(1)));
            % shift time course to start at zero
            t = t - t(1);
            intvl = [t(1) t(end)];
            cs = colormap(lines(64));
            
            subplot(421);

                plot(t,u); 
                title(['$(a)\qquad$Input signals $\{u^i\}_{i=1}^{' num2str(size(u,1)) '}$. $\Omega = 2\pi\cdot$' num2str(f) 'rad/s, $S = ' num2str(L) '\quad$']);
                xlim(intvl); ylabel('Amplitude');
                legend_text = cell(size(u,1),1);
                for i=1:size(u,1)
                    legend_text{i} = ['$u^' num2str(i) '(t)$'];
                end
                legend(legend_text,'location','east'); set(gca,'xticklabel',[]);

            subplot(423);
                plot(t,v(1,:)+bias,intvl, bias*ones(1,2),'--r');
                h_title = title(''); set(h_title,'Interpreter','latex'); title('$(b)\qquad$Biased filter output $v^1(t)+b$');
                hlabel = ylabel('Amplitude'); set(hlabel,'Interpreter', 'latex');
                h_legend = legend('$(u^1\ast h)(t)+b\quad$', ['Bias $b = ' num2str(bias, '%3.2f') '$'],'location','southeast'); set(h_legend,'Interpreter','latex');
                xlim([t(1) t(end)]); set(h_title,'Interpreter','latex'); set(gca,'xticklabel',[])

            subplot(425);
                plot(t,y1(1,:),t(logical(s(1,:))),y1(1,logical(s(1,:))),'.r', 'Markersize', 15);
                ylabel('$x^1$'); set(gca,'xticklabel',[])
                title('$(c)\qquad$ van der Pol oscillator output $x^1(t)$');
                xlim(intvl); ylim([-3 3]);
                legend(['$\mu=$' num2str(mu,'%3.1f')], 'Spikes (zero phase)$\quad$', 'location', 'east');

            subplot(427);
                for i=1:size(u,1)
                    stem(t(logical(s(i,:))) , 1.1*ones(1,sum(s(i,:))), '^', 'fill', 'color',cs(i,:) ); hold on;
                end
                box on; set(gca,'yticklabel',[],'ytick',[]);xlim(intvl);ylim([0 1.2]);
                title(['$(d)\qquad$Output sequences $(t^i_k)_{k=1}^{n^i}$. $n=\sum_{i=1}^{' num2str(size(u,1)) '}n^i = ' num2str(sum(sum(s)),'%3.0f') '\quad$']);
                legend_text = cell(size(u,1),1);
                for i=1:size(u,1)
                    legend_text{i} = ['$(t_k^' num2str(i) ')_{k=1}^{' num2str(sum(s(i,:))) '}$'];
                end
                legend(legend_text,'Location','East');xlabel('Time, [s]');xlim(intvl);ylim([0 1.2]);


            subplot(4,2,[2 4]);
                plot(lc(1,:), lc(2,:),'-m',y1(1,1:10:end),y2(1,1:10:end),...
                     lc(1,1), lc(2,1),'.r', 'Markersize', 15);
                xlabel('$x^1$');ylabel('$y^1$');
                title('$(e)\qquad$[Filter]-[van der Pol] response in the phase plane');
                legend(['Stable orbit, $b=' num2str(bias) '\quad$'],...
                       'Perturbed orbit$\quad$', 'Zero phase (spikes)');


            subplot(426);
                plot(t_filt, h, '--k', t_filt, Ph, '-b', t_filt, h_rec, '-r');
                axis([t_filt(1) t_filt(end) 1.1*min(h) 1.1*max(h)]);
                legend(['$h,\,$MSE$(\mathcal{P}h^*,h)$ =' num2str(10*log10(h_hhat_rmse^2), '%3.1f') 'dB $\quad$'],...
                       ['$\mathcal{P}h,\,$MSE$(\mathcal{P}h^*,\mathcal{P}h)$ =' num2str(10*log10(ph_hhat_rmse^2), '%3.1f') 'dB $\qquad$'],...
                       ['$\mathcal{P}h^*$, from $n = ' num2str(sum(sum(s))) '$ spikes'],...
                       'Location','NorthEast');
                title('$(e)\qquad$Original filter vs. the identified filter');
                ylabel('Amplitude');xlabel('Time, [s]');

            subplot(428);
                t_fft = dt*(-2*numel(t_filt):2*numel(t_filt));
                h_rec_rep = repmat(h_rec(1:end-1),1,4);
                n_fft = 2^(1+nextpow2(numel(h_rec_rep)));
                h_long = [ zeros(1,numel(find(t_fft<t_filt(1)))) h zeros(1,numel(find(t_fft>t_filt(end))))];
                [Fh freq_h] = cim_utility.centeredFFT(h_long, 1/dt, n_fft);
                [Fh_rec freq_h_rec] = cim_utility.centeredFFT(h_rec_rep, 1/dt, n_fft);
                plot(freq_h_rec, 10*log10(abs(Fh_rec)),'-b',...
                     freq_h, 10*log10(abs(Fh)),'-r' );
                xlim([-125 125]);ylim([-40 10]);grid on;
                ylabel('$10\log|\mathcal{F}\hat{h}|\;$');xlabel('Frequency, [Hz]');
                legend('supp$(\mathcal{FP}h^*)=[-\Omega,\,\Omega]\quad$',...
                       'supp$(\mathcal{F}h)\supset[-\Omega,\,\Omega]\quad$',...
                       'location','southeast');
                title('$(g)\qquad$Fourier amplitude spectra of $\mathcal{P}h^*$ and $h$');
        end
        
        function plotPaperFig9(t_filt,dt,h,Ph_cell,f_v,T)
            N = 2^(1+nextpow2(length(t_filt)));

            Fs = 1/dt;                                          % get the sampling frequency
            [Fh freq_h] = cim_utility.centeredFFT(h, Fs, N);

            for i=1:length(f_v)
                for j = 1:size(T,2)
                    L = f_v(i)*T(j);
                    l_idx = j + (i-1)*4;
                    r_idx = (j+2) + (i-1)*4;
                    [FPh freq_Ph] = cim_utility.centeredFFT(Ph_cell{i,j}, Fs, N);

                    subplot(3,4,l_idx);
                        plot(t_filt, h,'-r',t_filt, Ph_cell{i,j},'-b');
                        xlim([-0.3 0.3]);
                        legend('$h$', ['$\mathcal{P}h$, $L = ' num2str(L) '\;$']);
                    subplot(3,4, r_idx);
                        plot(freq_Ph,10*log10(abs(FPh)),'-b',...
                             freq_h, 10*log10(abs(Fh)),'-r');
                        xlim([-125 125]); grid on; ylim([-60 20]);
                        legend('$10\log|\mathcal{F}(h)|$',...
                               ['$10\log|\mathcal{F}(\mathcal{P}h)|$, $L = ' num2str(L) '\quad$'],...
                               'location','southeast');
                    if i==1
                        subplot(3,4,l_idx);
                            title(['$T = ' num2str(L/f_v(i)) '\,$s']);
                        subplot(3,4,r_idx);
                            title(['$T = ' num2str(L/f_v(i)) '\,$s']);
                    end
                    if j==1
                        subplot(3,4,l_idx);
                            ylabel(['$\Omega = 2\pi\cdot' num2str(f_v(i)) '\,$rad/s']);
                    else
                        subplot(3,4,l_idx);set(gca,'yticklabel',[]);
                        subplot(3,4,r_idx);set(gca,'yticklabel',[]);
                    end
                    if i==length(f_v)
                        subplot(3,4,l_idx);xlabel('Time, [s]');
                        subplot(3,4,r_idx);xlabel('Frequency, [Hz]');
                    else
                        subplot(3,4,l_idx);set(gca,'xticklabel',[]);
                        subplot(3,4,r_idx);set(gca,'xticklabel',[]);
                    end
                end
            end
        end
        function plotPaperFig11(T1,T2,f,L,bias,delta,dt,t_filt,h,Ph,t,u,v,...
                         z,voltage,h_rec,total_spikes_used)
            idx = find( t_filt >= T1 & t_filt <= T2 );
            h_hhat_err   = diag(1./max(abs(h),[],2))*(h-h_rec);
            ph_hhat_err  = diag(1./max(abs(Ph),[],2))*(Ph-h_rec);
            h_hhat_rmse  = sqrt( dt*trapz(h_hhat_err(:,idx).^2,2)/dt/numel(idx));
            ph_hhat_rmse = sqrt( dt*trapz(ph_hhat_err(:,idx).^2,2)/dt/numel(idx));
            intvl = [t(1) t(end)];

            cs = colormap(lines(64));
            subplot(421);
                plot( t, u(:,:,1)); hold on; 
                title(['$(a)\qquad$Input signal triplet $\mathbf{u}^1$. $\Omega = 2\pi\cdot$'...
                      num2str(f) 'rad/s, $L = ' num2str(L) '\quad$']);
                xlim(intvl); ylabel('Amplitude');set(gca,'xticklabel',[]);
                legend_text = cell(size(u,1),1);
                for i=1:size(u,1)
                    legend_text{i} = ['$u^{1' num2str(i) '}(t)$'];
                end
                legend(legend_text,'location','east');

            subplot(423);
                plot(t, v(1,:)-z(1,:), intvl, bias*ones(1,2),'--r', ...
                     intvl, zeros(1,2),'--k', intvl, -bias*ones(1,2),'--r');
                xlim(intvl); set(gca,'xticklabel',[]); ylabel('Amplitude'); 
                title('$(b)\qquad$Biased filter output $\sum_{m=1}^3(u^{1m}\ast h^m)(t)-z^1(t)$');
                legend('$v^1(t)-z^1(t)$', ['Biases $\pm b$, $b = ' num2str(bias, '%3.2f') '\quad$'], 'Zero line','location','East');

            subplot(425);
                s = [diff(z(1,:)) 0]/2/bias;
                idx = s ~= 0;
                plot(intvl,delta*ones(1,2),'--r', t,voltage(1,:), '-b',...
                     t(idx), delta*s(idx), '.r', 'markersize',15);
                xlim(intvl);ylim([-1.2*delta 1.2*delta]);ylabel('Amplitude'); set(gca,'xticklabel',[]);
                legend(['Thresholds $\pm\delta$, $\delta = ' num2str(delta,'%3.1e') '\qquad$'],...
                        '$\int_{t_k}^t[v^1(s)-z^1(s)]ds,\,\forall k\quad$','Trigger times', 'Location','East');
                hold on; plot(intvl,-delta*ones(1,2),'--r'); 
                title('(c)$\qquad$ASDM response to the triplet $\mathbf{u}^1$');

            subplot(427);
                plot(t, z(1,:),intvl, zeros(1,2),'--k');hold on;
                h_spk = stem(t(idx), 0.75*bias*ones(1,sum(idx)), '^k', 'filled'); 
                set(get(h_spk,'BaseLine'),'LineStyle','--')
                xlim(intvl); ylim([-1.2*bias 1.2*bias]); ylabel('Amplitude'); xlabel('Time, [s]');
                title(['(d)$\qquad$ASDM output $z^1(t)$ and its zero crossings $(t_k^1)_{k=1}^{' num2str(sum(idx)) '}$']);
                legend('$z^1(t)$', 'Zero line', ['Spikes $(t_k^1)_{k=1}^{' num2str(sum(idx)) '}\quad$'], 'location','East');

            subplot(422);
                legend_text = cell(size(u,1),1);
                for i=1:size(z,1)
                    s_idx = find( diff(z(i,:)) ~= 0 );
                    hold on; stem(t(s_idx), 1.1*ones(1,numel(s_idx)),'^', 'fill', 'color',cs(i,:) );
                    legend_text{i} = ['$(t_k^' num2str(i) ')_{k=1}^{' num2str(numel(s_idx)) '}\;$'];
                end
                set(gca,'yticklabel',[],'ytick',[],'box', 'on');
                legend(legend_text,'location','east');
                title(['$(e)\qquad$Output sequences $(t^i_k)_{k=1}^{n^i}$. $n=\sum_{i=1}^{' num2str(size(u,1)) '}n^i = ' num2str(total_spikes_used,'%3.0f') '\quad$']);
                xlim(intvl); ylim([0 1.2]); set(gca,'xticklabel',[]);

            subplot(4,2,4);
                plot(t_filt, h(1,:),'--k',t_filt, Ph(1,:),'k');
                hold on;plot(t_filt, h_rec(1,:),'Color',cs(3,:));
                legend(['$h^1,\,$MSE$(\mathcal{P}h^{1*},h^1)$ =' num2str(10*log10(h_hhat_rmse(1)^2), '%3.1f') 'dB $\qquad$'],...
                       ['$\mathcal{P}h^1,\,$MSE$(\mathcal{P}h^{1*},\mathcal{P}h^1)$ =' num2str(10*log10(ph_hhat_rmse(1)^2), '%3.1f') 'dB $\qquad$'],...
                       ['$\mathcal{P}h^{1*}$, from $n = ' num2str(total_spikes_used) '$ spikes'],'Location','NorthEast');
                title('$(f)\qquad$Original filter $h^1$ vs. the identified filter $\mathcal{P}h^{1*}$');
                ylabel('Amplitude');xlim([t_filt(1) t_filt(end)]); ylim([-50 100]); set(gca,'xticklabel',[]);

            subplot(4,2,6);
                plot(t_filt, h(2,:),'--k',t_filt, Ph(2,:),'k');
                hold on;plot(t_filt, h_rec(2,:),'Color',cs(2,:));
                legend(['$h^2,\,$MSE$(\mathcal{P}h^{2*},h^2)$ =' num2str(10*log10(h_hhat_rmse(2)^2), '%3.1f') 'dB $\qquad$'],...
                       ['$\mathcal{P}h^2,\,$MSE$(\mathcal{P}h^{2*},\mathcal{P}h^2)$ =' num2str(10*log10(ph_hhat_rmse(2)^2), '%3.1f') 'dB $\qquad$'],...
                       ['$\mathcal{P}h^{2*}$, from $n = ' num2str(total_spikes_used) '$ spikes'],'Location','NorthEast');
                title('$(g)\qquad$Original filter $h^2$ vs. the identified filter $\mathcal{P}h^{2*}$');
                ylabel('Amplitude');xlim([t_filt(1) t_filt(end)]); ylim([-50 100]); set(gca,'xticklabel',[]);

            subplot(4,2,8);
                plot(t_filt, h(3,:),'--k',t_filt, Ph(3,:),'k');
                hold on;plot(t_filt, h_rec(3,:),'Color',cs(1,:));
                legend(['$h^3,\,$MSE$(\mathcal{P}h^{3*},h^3)$ =' num2str(10*log10(h_hhat_rmse(3)^2), '%3.1f') 'dB $\qquad$'],...
                       ['$\mathcal{P}h^3,\,$MSE$(\mathcal{P}h^{3*},\mathcal{P}h^3)$ =' num2str(10*log10(ph_hhat_rmse(3)^2), '%3.1f') 'dB $\qquad$'],...
                       ['$\mathcal{P}h^{3*}$, from $n = ' num2str(total_spikes_used) '$ spikes'],'Location','SouthEast');
                title('$(h)\qquad$Original filter $h^3$ vs. the identified filter $\mathcal{P}h^{3*}$');
                ylabel('Amplitude');xlim([t_filt(1) t_filt(end)]); ylim([-100 50]);xlabel('Time, [s]');        
        end
        function plotPaperFig12(T1,T2,tau_1,tau_2,f,bias,delta,dt,t_filt,h,Ph,t,u,v,...
                         s,voltage,h_rec,total_spikes_used)
            idx = find( t_filt >= T1 & t_filt <= T2 );
            h_hhat_err   = diag(1./max(abs(h),[],2))*(h-h_rec);
            ph_hhat_err  = diag(1./max(abs(Ph),[],2))*(Ph-h_rec);
            h_hhat_rmse  = sqrt( dt*trapz(h_hhat_err(:,idx).^2,2)/dt/numel(idx));
            ph_hhat_rmse = sqrt( dt*trapz(ph_hhat_err(:,idx).^2,2)/dt/numel(idx));
            intvl = [tau_1 tau_2];

            cs = colormap(lines(64));
            subplot(421);
                plot( t, u); hold on; 
                title(['$(a)\qquad$Input signals $u^i\in\Xi$,$i=$1,...,' num2str(size(u,1))...
                       '. $\Omega = 2\pi\cdot$' num2str(f) 'rad/s$\quad$']);
                xlim(intvl); ylabel('Amplitude');set(gca,'xticklabel',[]);
                legend_text = cell(size(u,1),1);
                for i=1:size(u,1)
                    legend_text{i} = ['$u^{' num2str(i) '}$'];
                end
                legend(legend_text,'location','east');

            subplot(423);
                plot(t, v(1,:)+bias, intvl, bias*ones(1,2),'--r', ...
                     intvl, zeros(1,2),'--k');
                xlim(intvl); ylim([-1 1]);set(gca,'xticklabel',[]); ylabel('Amplitude'); 
                title('$(b)\qquad$Biased filter output $v^1(t)+b$');
                legend('$(u^1\ast h)(t)+b\quad$', ...
                       ['Bias $b = ' num2str(bias, '%3.2f') '$'], ...
                       'Zero line','location','SouthEast');

            subplot(425);
                plot(t, voltage(1,:),'-b',intvl, delta*ones(1,2),'--r');hold on;
                plot(t(logical(s(1,:))), delta*ones(1,sum(s(1,:))), 'or',...
                    'linewidth', 1, 'markersize',6, 'markerfacecolor','r');
                xlim(intvl);ylim([0 1.2*delta]);ylabel('Amplitude'); set(gca,'xticklabel',[]);
                legend('$\int_{t_k}^{t}(u^1\ast h)(s)ds,\,\forall k\qquad$',...
                       ['$\delta = ' num2str(delta) '$'], ...
                       '$\int_{t_k}^{t}(u^1\ast h)(s)ds=\delta$', 'Location','southeast');
                title('(c)$\qquad$Ideal IAF neuron response to $\mathbf{u}^1$');

            subplot(427);
                legend_text = cell(size(u,1),1);
                for i=1:size(u,1)
                    tk_idx = s(i,:)==1;
                    legend_text{i} = ['$(t_k^' num2str(i) ')_{k=1}^{' num2str(sum(tk_idx)) '}$'];
                    hold on;stem(t(tk_idx) , 1.1*ones(1,sum(tk_idx)),...
                                 '^', 'fill', 'color', cs(i,:) ); 
                end
                set(gca,'yticklabel',[],'ytick',[],'box','on');
                xlabel('Time, [s]');xlim(intvl);ylim([0 1.2]);
                title(['$(d)\qquad$Output sequences $(t^i_k)_{k=1}^{n^i}$. $\sum_{i=1}^{' num2str(size(u,1)) '}n^i = ' num2str(total_spikes_used,'%3.0f') '\quad$']);
                legend(legend_text,'Location','East');
                
            subplot(422);
                plot(t_filt, h, '--k', t_filt, Ph, '-b', t_filt, h_rec, '-r');
                axis([t_filt(1) t_filt(end) 1.1*min(h) 1.1*max(h)]);
                legend(['$h,\,$MSE$(\mathcal{P}h^*,h)$ =' num2str(10*log10(h_hhat_rmse^2), '%3.1f') 'dB $\quad$'],...
                       ['$\mathcal{P}h,\,$MSE$(\mathcal{P}h^*,\mathcal{P}h)$ =' num2str(10*log10(ph_hhat_rmse^2), '%3.1f') 'dB $\qquad$'],...
                       ['$\mathcal{P}h^*$, from $n = ' num2str(total_spikes_used) '$ spikes'],...
                       'Location','NorthEast');
                title('$(e)\qquad$Original filter vs. the identified filter');
                ylabel('Amplitude');xlabel('Time, [s]');

            subplot(4,2,4);
                t_fft = dt*(-2*numel(t_filt):2*numel(t_filt));
                K = 2*f*sinc(2*f*(t_fft));
                n_fft = 2^(1+nextpow2(numel(K)));
                [FK freq_K] = cim_utility.centeredFFT(K, 1/dt, n_fft);
                plot(freq_K, 10*log10(abs(FK)));
                title('$(f)\qquad$Fourier amplitude spectrum of $K$');grid on;
                legend('supp$(\mathcal{F}K)=[-\Omega,\,\Omega]\quad$','location','southeast');
                ylabel('$10\log|\mathcal{F}K|\;$'); set(gca,'xticklabel',[]);
                xlim([-125 125]);ylim([-40 20]);
            
            subplot(426);
                %h_pad = [ zeros(1,numel(find(t_fft<t_filt(1)))) h zeros(1,numel(find(t_fft>t_filt(end))))];
                [Fh freq_h] = cim_utility.centeredFFT(h, 1/dt, n_fft);
                plot(freq_h, 10*log10(abs(Fh)));
                legend('supp$(\mathcal{F}h)\supset[-\Omega,\,\Omega]\quad$','location','southeast');
                set(gca,'xticklabel',[]);grid on;
                xlim([-125 125]);ylim([-40 20]);ylabel('$10\log|\mathcal{F}h|\;$');
                title('$(g)\qquad$Fourier amplitude spectrum of $h$');

            subplot(428);
                [Fh_rec freq_h_rec] = cim_utility.centeredFFT(h_rec, 1/dt, n_fft);
                plot(freq_h_rec, 10*log10(abs(Fh_rec)));
                xlim([-125 125]);ylim([-40 20]);grid on;
                ylabel('$10\log|\mathcal{F}\mathcal{P}h^*|\;$');xlabel('Frequency, [Hz]');
                title('$(h)\qquad$Fourier amplitude spectrum of $\mathcal{P}h^*$');
                legend('supp$(\mathcal{F}\mathcal{P}h^*)=[-\Omega,\,\Omega]\quad$','location','southeast');
        
        end
    end
end
