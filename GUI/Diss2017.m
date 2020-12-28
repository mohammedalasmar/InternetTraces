function varargout = Diss2017(varargin)

% DISS2017 MATLAB code for Diss2017.fig
% Last Modified by GUIDE v2.5 29-Sep-2017 23:26:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Diss2017_OpeningFcn, ...
                   'gui_OutputFcn',  @Diss2017_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Diss2017 is made visible.
function Diss2017_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Diss2017 (see VARARGIN)

% Choose default command line output for Diss2017
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Diss2017 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Diss2017_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc    
% Get default command line output from handles structure
varargout{1} = handles.output;

set(handles.text62,'String',' ') ;
textstring=get(handles.text62,'string'); 
textstring=strcat(textstring,'Welcome');
set(handles.text62,'string',textstring);
 
axes(handles.axes5) , imshow('matlab_logo.jpg')
axes(handles.axes3) , imshow('2.jpg')
axes(handles.axes4) , imshow('11.jpg')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      R/S TEST       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global nn
global nnn

global ave VAR epsi T kk 
sto=0;y=0;z=1;R=0;

%x=[0.01 0.05 0.1 0.5 1 1.5 2 2.5 3.5 4 5 ];  % time scale value
x=str2num(get(handles.edit2,'string'));              % times of excution   

VAR=zeros(1,length(x)); C1=zeros(1,length(x)); C2=zeros(1,length(x));C3=C1;C4=C1;C5=C1;C6=C1;
ave=zeros(1,length(x)); 

for kk=1:length(x); 
sto=0;y=0;z=1;R=0;
T=x(kk);
for i=1:length(nn)                
    if nn(i)>=x(kk)
        sto(z)=i;
        x(kk)=x(kk)+T;
        z=z +1;
    end
end  
sto=[1 sto length(nn)+1];

for j=1:length(sto)-1  
    y(j)=sum(nnn(sto(j):sto(j+1)-1)); 
end
yy=y/T;  %% yy==>bit/sec

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segments=[ 2 4 8 16  32 64 128 256 512 1024 2048] ;
for ii=1:19;
    Seggg(ii)=2^(ii);
end

p=max(find(length(y)-Seggg>0));
Segments=Seggg(1:p-2);
ALL= Seggg(1:p);
for Seg=1:length(Segments); 
    NUMMM=ALL(end);
    DATAPOINTS=y(1:NUMMM);
    
    for w=1:Segments(Seg);
        AAA=DATAPOINTS(((NUMMM/Segments(Seg))*(w-1))+1:(NUMMM/Segments(Seg))*w)  ;
        AVERAGE(w)=mean(AAA);
        STDdev(w)=std(AAA);
        deviation= AAA-AVERAGE(w); 
        X=zeros(1,NUMMM/Segments(Seg));
        X(1)=deviation(1);
        
        for uu=2:length(X);
            X(uu)=X(uu-1)+deviation(uu);       % X: Cumulative sum of the deviation from the mean
        end        
        MAXvalue(w)=max(X);
        MINvalue(w)=min(X);
        Range(w)=MAXvalue(w)-MINvalue(w);
        R_S(w)=Range(w)/STDdev(w) ; 
        LOG_R_S(w)=log10(R_S(w));
       Variance(w)=var(AAA) ;
        
    end
    
    Xm=AVERAGE;
    VVVaarr(Seg)=var(Xm);
    Variance_tot(Seg)=sum(Variance) ;
    LOG_R_Stot(Seg)=mean(LOG_R_S);
    LOG_Size(Seg)=log10(NUMMM/Segments(Seg));
 
end
Variance_tot=Variance_tot/Variance_tot(length(Segments))

% H=1+0.5*slope   for time variance
RS_all_results=[Segments;10.^LOG_Size;LOG_Size;LOG_R_Stot];

Fitting= polyfit(LOG_Size,LOG_R_Stot,1); %% y=mx+b , to find the value of m and b 
Hurst_RS=Fitting(1);
% axes(handles.axes2) ,hold off
% for plotting
figure,
f_val=polyval(Fitting,LOG_Size);           %% the exact fitting line 
 plot(LOG_Size,LOG_R_Stot,'*',LOG_Size,f_val,'-r')

y1=LOG_Size+Fitting(2);     % y=1x+b
y2=0.5*LOG_Size+Fitting(2); % y=0.5x+b
hold on ,plot(LOG_Size,y1,'k','LineWidth', 2),hold on ,plot(LOG_Size,y2,'k','LineWidth', 2)
% Labels 
text(3.5, 3.3 , 'H=1', 'Color', 'k'),text(3.5, 1.2, 'H=0.5', 'Color', 'k');
text(1.2, 2.8 , num2str(Hurst_RS), 'Color', 'k')
xlabel('log(n)'),ylabel('log(R/S)'),title('Rescaled-Range');

set(handles.text4,'string',Hurst_RS); 




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      T-V  TEST       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
global y DATAPOINTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segments=[ 2 4 8 16  32 64 128 256 512 1024 2048] ;
for ii=1:19;
    Seggg(ii)=2^(ii);
end

p=max(find(length(y)-Seggg>0));
Segments=Seggg(1:p-2);
ALL= Seggg(1:p);
for Seg=1:length(Segments); 
    NUMMM=ALL(end);
    DATAPOINTS=y(1:NUMMM);
    
    for w=1:Segments(Seg);
        AAA=DATAPOINTS(((NUMMM/Segments(Seg))*(w-1))+1:(NUMMM/Segments(Seg))*w)  ;
        AVERAGE(w)=mean(AAA);
        STDdev(w)=std(AAA);
        deviation= AAA-AVERAGE(w); 
        X=zeros(1,NUMMM/Segments(Seg));
        X(1)=deviation(1);
        
        for uu=2:length(X);
            X(uu)=X(uu-1)+deviation(uu);       % X: Cumulative sum of the deviation from the mean
        end        
        MAXvalue(w)=max(X);
        MINvalue(w)=min(X);
        Range(w)=MAXvalue(w)-MINvalue(w);
        R_S(w)=Range(w)/STDdev(w) ; 
        LOG_R_S(w)=log10(R_S(w));
       Variance(w)=var(AAA) ;
    end
    
    Xm=AVERAGE;
    VVVaarr(Seg)=var(Xm);
    Variance_tot(Seg)=sum(Variance) ;
    LOG_R_Stot(Seg)=mean(LOG_R_S);
    LOG_Size(Seg)=log10(NUMMM/Segments(Seg));
 
end
Variance_tot=Variance_tot/Variance_tot(length(Segments))

% H=1+0.5*slope   for time variance

RS_all_results=[Segments;10.^LOG_Size;LOG_Size;LOG_R_Stot];
Fitting= polyfit(LOG_Size,LOG_R_Stot,1); %% y=mx+b , to find the value of m and b 
Hurst_RS=Fitting(1);
%%%% time variance test plotting and fitting
% Variance_tot=fliplr(Variance_tot)
Fitting_var_test= polyfit(LOG_Size,log10(VVVaarr),1);
Hurst_VARtest=1+0.5*Fitting_var_test(1);
f_val_timeVARtest=polyval(Fitting_var_test,LOG_Size); 

% axes(handles.axes2) ,hold off
figure,
plot(LOG_Size,log10(VVVaarr),'*',LOG_Size,f_val_timeVARtest,'-r')
y1_2=-1*LOG_Size+Fitting_var_test(2);     % y=1x+b
hold on ,plot(LOG_Size,y1_2,'k','LineWidth', 2)
text(3.5, 3.3 , 'Beta=-1', 'Color', 'k')
text(1.2, 2.8 , num2str(Hurst_RS), 'Color', 'k')
xlabel('log(m)'),ylabel('log(Variance)'),title('Variance-Time');

set(handles.text6,'string',Hurst_VARtest); 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Generate   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
global y  yy
global nn
global nnn
global epsi T

y=0;

x=str2num(get(handles.edit2,'string'));              % times of excution   
epsi=str2num(get(handles.edit5,'string'));              % 



for kk=1:length(x); 
    sto=0;y=0;z=1;
    T=x(kk);
    
    for i=1:length(nn)                
        if nn(i)>=x(kk)
            sto(z)=i;
            x(kk)=x(kk)+T;
            z=z +1;
        end
    end
    sto=[1 sto length(nn)+1];
    
    for j=1:length(sto)-1  
        y(j)=sum(nnn(sto(j):sto(j+1)-1)); 
    end
    yy=y/T;  %% yy==>bit/sec
end


disp('done')
figure(1)
c=0:T:nn(end);
subplot(length(x),1,kk), stairs(c,y/T) , xlabel('time (sec)'), ylabel('data rate (bps)'),title(num2str(T)) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Peridogram  TEST       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton6_Callback(hObject, eventdata, handles)

global   DATAPOINTS
%  sequence=nnn';

sequence=DATAPOINTS;

N= length(sequence);
Xk = fft(sequence) ;
I = abs(Xk).^2/(2*pi*N);

k=1:N;
w=2*pi*k/N;
 
% Use the lowest 20% part of periodogram to estimate the similarity.
W = w(1:floor(N*0.2));
II = I(1:floor(N*0.2));

WLog=log10(W);   IILog=log10(II);
 p1 = polyfit(WLog,IILog,1);
Yfit = polyval(p1,WLog); 
length(WLog)

length(Yfit)

HURST_perido=(1-p1(1))/2;

% axes(handles.axes2) ,hold off
figure,
plot(log10(w),log10(I),'b.'),hold on , plot(WLog,Yfit,'r')
xlabel('log10(Frequency)'),ylabel('log10(Periodogram)'),title('Periodogram Method');
 
set(handles.text8,'string',HURST_perido); 



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function text10_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to text10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      QQ plot      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
global y
global T

figure,
subplot(211),gqqplot(y(1:end-1)/T,'normal') ,grid on
%text(-1,-1 , num2str(gamma1(1,2)), 'Color', 'k')
subplot(212),histfit(y(1:end-1)/T,200,'normal'), 
xlabel('A(T)  bit/sec' ,'FontSize', 12), 
ylabel('PDF of A(T)', 'FontSize', 12),grid on
legend('A(T) distr.','Normal distr.')

y2=sort(y(1:end-1));   %% the sorted y
reff_Gauss=sort(randn(1,length(y)-1));
gamma_Gauss=corrcoef(y2,reff_Gauss);
gamma_Gauss1=gamma_Gauss(1,2);
set(handles.text29,'string',gamma_Gauss1); 



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
global y
global T

figure,
subplot(211),gqqplot(y(1:end-1)/T,'Lognormal')  ,grid on
subplot(212),histfit(y(1:end-1)/T,200,'Lognormal'),grid on
xlabel('A(T)  bit/sec' ,'FontSize', 12), 
ylabel('PDF of A(T)', 'FontSize', 12),
legend('A(T) distr.','Lognormal distr.')

%%% Correlation coefficient:: Lognormal
y2=sort(y(1:end-1));   %% the sorted y
reff_Lognorm=sort(lognrnd(0,1,[1 length(y)-1]));
gamma_Lognorm=corrcoef(y2,reff_Lognorm);
gamma_Lognorm1=gamma_Lognorm(1,2);
set(handles.text25,'string',gamma_Lognorm1); 


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
global y
global T
figure,
subplot(211),gqqplot(y(1:end-1)/T,'gev')  ,grid on
subplot(212),histfit(y(1:end-1)/T,200,'gev'),grid on

xlabel('A(T)  bit/sec' ,'FontSize', 12), 
ylabel('PDF of A(T)', 'FontSize', 12),
legend('A(T) distr.','GEV distr.')

%%% Correlation coefficient::Generalized Extreme Value
y2=sort(y(1:end-1));   %% the sorted y
reff_GEV=sort(gevrnd(0,1,0,[1 length(y)-1]));
gamma_GEV=corrcoef(y2,reff_GEV);
gamma_GEV1=gamma_GEV(1,2);
set(handles.text27,'string',gamma_GEV1); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Capacity     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)

global ave VAR   C1 C2 C3 C4 C5   

global T
global epsi
global y

%  STATISTICS Values 
VAR=var(y(1:end-1))  ; % bit^2
ave=mean(y(1:end-1))/T ; %% bit/sec
sigma=std(y(1:end-1));  %% bit

NN=VAR/T^2;
NN_root=sqrt(VAR/T^2);


%  Finding the Value of the Capacity
%%% first equ --- direct Gaussian way
% C1(kk)=2.332*sigma(kk)/T+ave(kk);
C1=icdf('norm',1-epsi,ave,sigma/T);

%%% second equ --- Nick equ
  C2=fsolve('solveC2',1.4*C1);

%%% third  equ --- C3-Meent et al. equ
C3=ave+(1/T)*sqrt(-2*log(epsi)*VAR);
% safety_margin=(1/T)*sqrt(-2*log(epsi)*VAR);
% safety_margin_part=(1/T)*sqrt(VAR);

%%% fourth equ --- Log-normal
% [mu,sigma]=lognfit(yy);
LOGnormalfit=lognfit(y/T);
C4=icdf('Lognormal',1-epsi,LOGnormalfit(1),LOGnormalfit(2));

%%% Fifth  equ --- Generalized Extreme Value
% [k_shape parameter,sigma,mu]=gevfit(yy);
GEVfit=gevfit(y/T);
C5=icdf('Generalized Extreme Value',1-epsi,GEVfit(1),GEVfit(2),GEVfit(3));

set(handles.text15,'string',C1); 
set(handles.text17,'string',C2); 
set(handles.text19,'string',C3); 
set(handles.text21,'string',C4); 
set(handles.text23,'string',C5); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Empirical C     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton11_Callback(hObject, eventdata, handles)
global    C1 C2 C3 C4 C5 
global   yy

epsi_c1=length(find(yy(1:end-1)>C1))/(length(yy)-1);   % direct 
epsi_c2=length(find(yy(1:end-1)>C2))/(length(yy)-1);    % nick
epsi_c3=length(find(yy(1:end-1)>C3))/(length(yy)-1);  % C3-Meent et al.
epsi_c4=length(find(yy(1:end-1)>C4))/(length(yy)-1); % lognormal
epsi_c5=length(find(yy(1:end-1)>C5))/(length(yy)-1); % GEV 


set(handles.text41,'string',epsi_c1); 
set(handles.text42,'string',epsi_c2); 
set(handles.text43,'string',epsi_c3); 
set(handles.text44,'string',epsi_c4); 
set(handles.text45,'string',epsi_c5); 


% --- Executes on button press in pushbutton12.



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
 close all
clc
clear all
% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
axes(handles.axes8) , imshow('report.jpg')
open('Report.pdf')


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
axes(handles.axes8)  , imshow('by.jpg')



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
global ave VAR   

global T

global y

%  STATISTICS Values 
VAR=var(y(1:end-1))  ; % bit^2
ave=mean(y(1:end-1))/T ; %% bit/sec
sigma=std(y(1:end-1));  %% bit

set(handles.text57,'string',ave); 
set(handles.text59,'string',VAR); 
set(handles.text61,'string',sigma); 


% --- Executes on button press in pushbutton19.



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
 global nn nnn 
set(handles.text62,'String',' ') ;
textstring=get(handles.text62,'string'); 
textstring=strcat(textstring,'Please wait ..');
set(handles.text62,'string',textstring);   

[filename1,filepath1]=uigetfile({'*.txt','All Files'},'Select Data File');
fullpathname=strcat(filepath1,filename1);
fid = fopen(fullpathname);


   
a = fscanf(fid,'%g %g',[2 inf]); % It has two rows now.


a=a' ;
nn=a(:,1);  % time vector
nnn=a(:,2)*8 ; % packets length 
    
set(handles.text62,'String',' ') ;
textstring=get(handles.text62,'string'); 
textstring=strcat(textstring,'You can start now ..');
set(handles.text62,'string',textstring);   
   
   
    
set(handles.text64,'string',length(nnn)); 
set(handles.text66,'string',sum(nnn)); 


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
close all

Diss2017
 


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
clc
clear all


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
Diss2017part2

function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 

 

 
