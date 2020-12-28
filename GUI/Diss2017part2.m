function varargout = Diss2017part2(varargin)
% DISS2017PART2 MATLAB code for Diss2017part2.fig
% Last Modified by GUIDE v2.5 29-Sep-2017 20:23:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Diss2017part2_OpeningFcn, ...
                   'gui_OutputFcn',  @Diss2017part2_OutputFcn, ...
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


% --- Executes just before Diss2017part2 is made visible.
function Diss2017part2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Diss2017part2 (see VARARGIN)

% Choose default command line output for Diss2017part2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Diss2017part2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Diss2017part2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
   
global nn nnn 
global ave VAR  T kk epsi  epsi_c1 epsi_c2 epsi_c3 epsi_c4 epsi_c5 T_array EPSI00
global C1  C2 C3 C4 C5 
T_array=[];
epsi_c1=[];
epsi_c2=[];
epsi_c3=[];
epsi_c4=[];
epsi_c5=[];
 
sto=0;y=0;z=1;R=0;

%     x=[0.01 0.05 0.1 0.5 1 1.5 2 2.5 3.5 4 5 ];  % time scale value
x=str2num(get(handles.edit1,'string'));              % times of excution   

EPSI00=str2num(get(handles.edit4,'string'));              % times of excution   

VAR=zeros(1,length(x)); C1=zeros(1,length(x)); C2=zeros(1,length(x));C3=C1;C4=C1;C5=C1;C6=C1;
ave=zeros(1,length(x)); epsi=zeros(1,length(x));

for kk=1:length(x); 
sto=0;y=0;z=1;R=0;
T=x(kk);
T_array(kk)=T ;

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
%  y=y(1:end-1);
%  finding the value of C1 and C2

yy=y/T;  %% yy==>bit/sec


%% STATISTICS Values 
VAR(kk)=var(y(1:end-1))  ; % bit^2
epsi(kk)=EPSI00  ;
ave(kk)=mean(y(1:end-1))/T ; %% bit/sec
sigma(kk)=std(y(1:end-1));  %% bit

NN(kk)=VAR(kk)/T^2;
NN_root(kk)=sqrt(VAR(kk)/T^2);
%% Finding the Value of the Capacity
%%% first equ --- C3-Meent et al.
C1(kk)=ave(kk)+(1/T)*sqrt(-2*log(epsi(kk))*VAR(kk));
safety_margin(kk)=(1/T)*sqrt(-2*log(epsi(kk))*VAR(kk));
safety_margin_part(kk)=(1/T)*sqrt(VAR(kk));


%%% second equ --- Nick equ
C2(kk)=fsolve('solveC2K',1.4*C1(kk));

%%% third equ --- direct Gaussian way
% C3(kk)=2.332*sigma(kk)/T+ave(kk);
C3(kk)=icdf('norm',1-epsi(kk),ave(kk),sigma(kk)/T);

%%% fourth equ --- Generalized Extreme Value
% [k_shape parameter,sigma,mu]=gevfit(yy);

GEVfit=gevfit(y/T);
C4(kk)=icdf('Generalized Extreme Value',1-epsi(kk),GEVfit(1),GEVfit(2),GEVfit(3));
% sum(pdf('Generalized Extreme Value',C4:1:2.6e7,GEVfit(1),GEVfit(2),GEVfit(3)))
% ans =0.0100

%%% Fifth equ --- Log-normal
% [mu,sigma]=lognfit(yy);

LOGnormalfit=lognfit(y/T);
   C5(kk)=icdf('Lognormal',1-epsi(kk),LOGnormalfit(1),LOGnormalfit(2));
%    C5(kk)=icdf('Lognormal',1-epsi(kk),ave(kk),sigma(kk)/T);

%% finding the correlation coefficient :: Gaussian Dis.
y2=sort(y(1:end-1));   %% the sorted y

reff_Gauss=sort(randn(1,length(y)-1));
gamma_Gauss=corrcoef(y2,reff_Gauss);
gamma_Gauss1(kk)=gamma_Gauss(1,2);


%%% Correlation coefficient::Generalized Extreme Value
reff_GEV=sort(gevrnd(0,1,0,[1 length(y)-1]));
gamma_GEV=corrcoef(y2,reff_GEV);
gamma_GEV1(kk)=gamma_GEV(1,2);

%%% Correlation coefficient:: Lognormal
reff_Lognorm=sort(gevrnd(0,1,0,[1 length(y)-1]));
gamma_Lognorm=corrcoef(y2,reff_Lognorm);
gamma_Lognorm1(kk)=gamma_Lognorm(1,2);

 
%% Empirical Epsi
epsi_c1(kk)=length(find(yy(1:end-1)>C1(kk)))/(length(yy)-1);   % C3-Meent et al.
epsi_c2(kk)=length(find(yy(1:end-1)>C2(kk)))/(length(yy)-1);    % nick
epsi_c3(kk)=length(find(yy(1:end-1)>C3(kk)))/(length(yy)-1);  % direct
epsi_c4(kk)=length(find(yy(1:end-1)>C4(kk)))/(length(yy)-1); % GEV
epsi_c5(kk)=length(find(yy(1:end-1)>C5(kk)))/(length(yy)-1); % lognormal

Corr_Coeff_All(kk,:)=[T,gamma_Gauss1(kk),gamma_GEV1(kk),gamma_Lognorm1(kk)]';
vvv(kk,:)=[T,epsi(kk),C3(kk),C2(kk),C1(kk),C5(kk),C4(kk),ave(kk),VAR(kk)/(T^2) ,sigma(kk)]';
Epsi_all(kk,:)=[T,epsi(kk),epsi_c3(kk),epsi_c2(kk),epsi_c1(kk),epsi_c5(kk),epsi_c4(kk)]';
end
  
set(handles.uitable1,'Data',vvv);
set(handles.uitable3,'Data',Epsi_all);

     
% --------------------------------------------------------------------
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
 global nn
global nnn
global C1  C2 C3 C4 C5 

global T kk 
sto=0;y=0;z=1;R=0;

% x=[1];
x=str2num(get(handles.edit1,'string'));              % times of excution   


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
    
    %%% Plotting 
    figure(1)
     c=0:T:nn(end);
    subplot(length(x),1,kk), 
    stairs(c,y/T) , xlabel('time (sec)'), ylabel('data rate (bps)'),title(['T= ', num2str(T)]) 
    hline3 = refline([0 C1(kk)]);
    hline3.Color = 'r';
  hline4 = refline([0 C5(kk)]);
    hline4.Color = 'g';
      hline5 = refline([0 C4(kk)]);
    hline5.Color = 'black';
 
legend([hline3 hline4 hline5],{'C3 ','C4 ','C5 '})
end
 

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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
  global nn nnn 
set(handles.text1,'String',' ') ;
textstring=get(handles.text1,'string'); 
textstring=strcat(textstring,'Please wait ..');
set(handles.text1,'string',textstring);   

[filename1,filepath1]=uigetfile({'*.txt','All Files'},'Select Data File');
fullpathname=strcat(filepath1,filename1);
fid = fopen(fullpathname);

   
a = fscanf(fid,'%g %g',[2 inf]); % It has two rows now.


a=a' ;
nn=a(:,1);  % time vector
nnn=a(:,2)*8 ; % packets length 
set(handles.text1,'String',' ') ;
textstring=get(handles.text1,'string'); 
textstring=strcat(textstring,'You can start now ..');
set(handles.text1,'string',textstring);   
   
   

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


% --- Executes on button press in pushbutton4.
 

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all



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


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global   epsi_c1 T_array EPSI00
err =epsi_c1;
y =EPSI00* ones(1,length(epsi_c1)); 

figure 
hold on
for i = 1:length(err)
    h=bar(i,err(i));
    if err(i) < EPSI00
        set(h,'FaceColor','g');
    else
        set(h,'FaceColor','r');
    end
end
hold off

% bar(err,'r')
 
for i=1:length(T_array)
text(i, err(i)-0.0004 ,num2str(err(i)));
end
  xticks([])
for i=1:length(T_array)
text(i-0.4, 0.001 , 'T=');
text(i-0.15 , 0.001 ,num2str(T_array(i))  , 'Color','blue');
end
grid  , hold on 
plot(y,'b--')
set(gca,'FontSize',14)

xlabel(' ','FontSize', 14)
ylabel('Emperical Epsi','FontSize', 16) 



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global   epsi_c2 T_array EPSI00
err =epsi_c2;
y =EPSI00* ones(1,length(epsi_c2)); 
figure 

figure 
hold on
for i = 1:length(err)
    h=bar(i,err(i));
    if err(i) < EPSI00
        set(h,'FaceColor','g');
    else
        set(h,'FaceColor','r');
    end
end
hold off
 
for i=1:length(T_array)
text(i, err(i)-0.0004 ,num2str(err(i)));
end
  xticks([])
for i=1:length(T_array)
text(i-0.4, 0.001 , 'T=');
text(i-0.15 , 0.001 ,num2str(T_array(i))  , 'Color','blue');
end
grid  , hold on 
plot(y,'b--')
set(gca,'FontSize',14)

xlabel(' ','FontSize', 14)
ylabel('Emperical Epsi','FontSize', 16) 


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global   epsi_c3 T_array EPSI00
err =epsi_c3;
y =EPSI00* ones(1,length(epsi_c3)); 
figure 
figure 
hold on
for i = 1:length(err)
    h=bar(i,err(i));
    if err(i) < EPSI00
        set(h,'FaceColor','g');
    else
        set(h,'FaceColor','r');
    end
end
hold off
 
for i=1:length(T_array)
text(i, err(i)-0.0004 ,num2str(err(i)));
end
  xticks([])
for i=1:length(T_array)
text(i-0.4, 0.001 , 'T=');
text(i-0.15 , 0.001 ,num2str(T_array(i))  , 'Color','blue');
end
grid  , hold on 
plot(y,'b--')
set(gca,'FontSize',14)

xlabel(' ','FontSize', 14)
ylabel('Emperical Epsi','FontSize', 16) 

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global   epsi_c5 T_array EPSI00
err =epsi_c5;
y =EPSI00* ones(1,length(epsi_c5)); 
figure 
figure 
hold on
for i = 1:length(err)
    h=bar(i,err(i));
    if err(i) < EPSI00
        set(h,'FaceColor','g');
    else
        set(h,'FaceColor','r');
    end
end
hold off
 
for i=1:length(T_array)
text(i, err(i)-0.0004 ,num2str(err(i)));
end
  xticks([])
for i=1:length(T_array)
text(i-0.4, 0.001 , 'T=');
text(i-0.15 , 0.001 ,num2str(T_array(i))  , 'Color','blue');
end
grid  , hold on 
plot(y,'b--')
set(gca,'FontSize',14)

xlabel(' ','FontSize', 14)
ylabel('Emperical Epsi','FontSize', 16) 


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global   epsi_c4 T_array EPSI00
err =epsi_c4;
y =EPSI00* ones(1,length(epsi_c4)); 
figure 
figure 
hold on
for i = 1:length(err)
    h=bar(i,err(i));
    if err(i) < EPSI00
        set(h,'FaceColor','g');
    else
        set(h,'FaceColor','r');
    end
end
hold off
 
for i=1:length(T_array)
text(i, err(i)-0.0004 ,num2str(err(i)));
end
  xticks([])
for i=1:length(T_array)
text(i-0.4, 0.001 , 'T=');
text(i-0.15 , 0.001 ,num2str(T_array(i))  , 'Color','blue');
end
grid  , hold on 
plot(y,'b--')
set(gca,'FontSize',14)

xlabel(' ','FontSize', 14)
ylabel('Emperical Epsi','FontSize', 16) 
