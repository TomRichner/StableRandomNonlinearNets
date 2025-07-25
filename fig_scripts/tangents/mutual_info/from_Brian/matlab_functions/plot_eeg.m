function plot_eeg(data,Fs,interval,varargin)

%plot_eeg(data,Fs,interval,varargin)
%data - rows-time points; col=channels
%interval - seconds per screen
%varargin - array of detected spikes for each channel

%INPUT: 
%Matrix of EEG data, where each channel is a column, negative is down
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultLineLineWidth',1)
set(0,'DefaultAxesFontSize',10)
set(0,'defaultAxesBox','off');
set(0,'defaultLineMarkerSize',12)
set(0, 'DefaultFigureColor', [.95 .95 .95]);
set(0,'defaultAxesFontSize',18)

global N scale t d offset t1 t2 f int h DetSpikes

%Plot EEG data
[M N] = size(data);
offset = 1:N;
f = Fs;
t = 1/Fs:1/Fs:M/Fs;
scale = 1;
d = data;
int = interval;
t1 = 1/f;
t2 = int;
if ~isempty(varargin)
    DetSpikes = varargin{1};
else
    DetSpikes = [];
end
h = figure;
draw_eeg()
xlim([t1 t2])


set(h,'KeyPressFcn',@KeyPressCb);

end

function s = KeyPressCb(~,evnt)
global N scale t d offset t1 t2 f int
%fprintf('key pressed: %s\n',evnt.Key);

if strcmp(evnt.Key,'uparrow')==1
    scale = scale*2;
    draw_eeg()
    xlim([t1 t2])
elseif strcmp(evnt.Key, 'downarrow')==1
    scale = scale/2;
    draw_eeg()
    xlim([t1 t2])
%elseif strcmp(evnt.Key,'space')==1

elseif strcmp(evnt.Key, 'rightarrow')==1
    t1=t2;
    t2=t2+int;
    draw_eeg()
elseif strcmp(evnt.Key, 'leftarrow')==1
    if t1>=int
        t2=t1;
        t1=(t1-int)+1/f;
        draw_eeg()
    end
end
end

function draw_eeg()
global N scale t d offset t1 t2 f h DetSpikes
figure(h)
clf

i1 = round(t1*f);
i2 = round(t2*f-1);
for i = 1:N %N is number of channels
    plot_data = d(i1:i2,i);
    plot_data = plot_data./500;
    plot_data = plot_data-mean(plot_data);
    plot(t(i1:i2),-scale*plot_data+offset(i)),hold on
    if ~isempty(DetSpikes)
        ySpikes = offset(i)*ones(1,length(DetSpikes{i}));
        plot(t(DetSpikes{i}),ySpikes,'k.')
    end
end
ylabel('Channels')
xlabel('Time')
set(gca,'Ytick',1:2:N)
ylim([0 N+1])
xlim([t1 t2])
grid on
end


