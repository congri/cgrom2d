%performance tests

clear Tf;

%load first Tf
tic
Tf = Tffile.Tf(:,1);
toc
clear Tf;

%load last Tf
tic
Tf = Tffile.Tf(:,1024);
toc
clear Tf;

%load 10 Tf's
tic
Tf = Tffile.Tf(:, 134:144);
toc
clear Tf;

%load all Tf's
tic
Tf = Tffile.Tf;
toc
clear Tf;