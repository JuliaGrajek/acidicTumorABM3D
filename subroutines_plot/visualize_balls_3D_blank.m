% JN Kather, NCT Heidelberg, 2017
% updated by J Grajek, 2022
function visualize_balls_3D_blank(mySystem,cnst,allSummaries, r)

currentSummary = allSummaries{end}; % current summary is last state
mysqueeze = @(varargin) (varargin); % squeeze function

for i=1:numel(allSummaries)
    plotTimeline(i).TU_Num = allSummaries{i}.TU_Num;
    plotTimeline(i).minpH = allSummaries{i}.minpH;
    plotTimeline(i).minO2 = allSummaries{i}.minO2;
    plotTimeline(i).minG = allSummaries{i}.minG;
    plotTimeline(i).minATP = allSummaries{i}.minATP;
    plotTimeline(i).PDL = allSummaries{i}.PDLfrac;
    plotTimeLine(i).CA9freq = allSummaries{i}.CA9freq;
    plotTimeLine(i).PDL1SuppProb = allSummaries{i}.PDL1SuppProb;
end

if numel(mySystem.TU.TUcells>0) %&& numel(mySystem.IM.IMcells>0)
    
    % prepare tumor antigenicity grid before cutting out tumor cell slice
    antigenicityGrid = double(zeros(size(mySystem.grid.Ln)));   % preallocate antigenicity grid
    antigenicityGrid(mySystem.TU.TUcells) = mySystem.TU.TUprop.Antigen;   % preallocate antigenicity grid
    
    %
    figure(1)    % prepare figure 1 = tumor 3D plot
    clf('reset')
    set(gcf, 'Position',  [100, 100, 200, 200])
    %     show full domain and add decorations
    xsize = 4; ysize = 10;
    posi.MainPanel = [1:6,11:16,21:26,31:36];
    posi.TimelinePanel = 7:8;
    posi.minpHPanel = 17:19;
    posi.spatialpHPanel = [27:30,37:40];
    posi.PDLPanel = 9:10;
    
    subplot(xsize,ysize,posi.MainPanel);
    
    TUcolors = double(hot(round(2*cnst.tumorColorLevels))); %colormap for tumor cells
    TUcolors = TUcolors(0.5*cnst.tumorColorLevels:1.5*cnst.tumorColorLevels,:);% crop colormap
    IMcolors = flipud(double(blugr(double(mySystem.params.IMkmax)+3))); % color map for lymphocytes
    
    N1 = mySystem.grid.N;
    N2 = mySystem.grid.M;
    N3 = mySystem.grid.P;
    
    %clearing cells from the boundary
    L = false(size(mySystem.grid.Ln));
    Lymph = L;
    Nec = L;
    Fib = L;
    
    clear90 = true;
    
    if clear90
        %clearing cells within 90 degree angle to have cutout
        [x, y, ~] = ind2sub([N1 N2 N3],mySystem.TU.TUcells);
        
        %calculating center of mass
        X = mean(x); Y = mean(y); %Z = mean(z);
        indx = (x > X) & (y > Y);
        mySystem.TU.TUcells(indx) = []; % remove tumor cells
        mySystem.TU.TUprop.Antigen(indx) = []; % also remove antigenicity property
        
        [x, y, ~] = ind2sub([N1 N2 N3],mySystem.IM.IMcells);
        mySystem.IM.IMcells((x > X) & (y > Y)) = [];
        
        necro = int32(find(mySystem.grid.Ln)');
        [x, y, ~] = ind2sub([N1 N2 N3],necro);
        necro((x > X) & (y > Y)) = [];
        %
        fibr = int32(find(mySystem.grid.Lf)');
        [x, y, ~] = ind2sub([N1 N2 N3],fibr);
        fibr((x > X) & (y > Y)) = [];
    end
    
    % update grids
    L( mySystem.TU.TUcells) = true;
    Lymph( mySystem.IM.IMcells) = true;
    
    Nec(necro) = true;
    Fib(fibr) = true;
    
    
    blLevels = 1;
    if blLevels %if perform blur
        
        %auxilary variable with indices to the cell neighborhood
        aux = int32([[-N1-1 -N1 -N1+1 -1 1 N1-1 N1 N1+1] ...
            [-N1-1 -N1 -N1+1 -1 1 N1-1 N1 N1+1]-N1*N2 ...
            [-N1-1 -N1 -N1+1 -1 1 N1-1 N1 N1+1]+N1*N2])';
        
        S = [mySystem.TU.TUcells unique(reshape(bsxfun(@plus,mySystem.TU.TUcells,aux),1,[]))];
        S2 = bsxfun(@plus,S,aux);
        S2(S2<1) = []; S2(S2>N1*N2*N3) = [];
        
        SN = [necro unique(reshape(bsxfun(@plus, necro,aux),1,[]))];
        SN2 = bsxfun(@plus,SN,aux);
        SN2(SN2<1) = []; SN2(SN2>N1*N2*N3) = [];
        
        %changing lattice from logical variable to float
        L = single(L);
        Lymph = single(Lymph);
        Nec = single(Nec);
        Fib = single(Fib);
        for i = 1:blLevels %for number of blurs
            L(S) = mean(L(S2)); %taking the average of neighborhood
            Nec(SN) = mean(Nec(SN2)); %taking the average of neighborhood
        end
    end
    
    %calculating isosurfaces and plotting
    hold on
    
    %necrosis
    pNf = patch(isosurface(1:N1,1:N2,1:N3,Nec,0.5));
    isonormals(1:N1,1:N2,1:N3,Nec,pNf)
    set(pNf,'FaceColor','k','EdgeColor','none');
    
    
    %fibrosis
    pLf = patch(isosurface(1:N1,1:N2,1:N3,Fib,0.5));
    isonormals(1:N1,1:N2,1:N3,Fib,pLf)
    set(pLf,'FaceColor',[128 128 128]/255,'EdgeColor','none');
    
    % lymphocytes
    fvLym = isosurface(1:N1,1:N2,1:N3,Lymph,0.25);
    pL = patch(fvLym);
    isonormals(1:N1,1:N2,1:N3,Lymph,pL)
    set(pL,'FaceColor',IMcolors(end,:),'EdgeColor','none');
    
    
    % plot tumor as one color only
    fvTum = isosurface(1:N1,1:N2,1:N3,L,0.25);
    pT = patch(fvTum);
    isonormals(1:N1,1:N2,1:N3,L,pT)
    isonormals(1:N1,1:N2,1:N3,L,pT)
    set(pT,'FaceColor',TUcolors(1,:),'EdgeColor','none');
    
    hold off
    
    set(gcf,'Color','w');
    %axis equal
    axis off
    
    xlim([1 N1]);
    ylim([1 N2]);
    zlim([1 N3]);
    
    view([66 43]); %135 20
    camlight
    camlight('right') % add more light
    camlight('left')  % add more light
    lighting gouraud
    
    
    %     prepare decorations
    myPos = get(0,'Screensize');
    myPos(1) = myPos(1) + myPos(3)*(1/10);
    myPos(2) = myPos(2) + myPos(4)*(1/10);
    myPos(3) = myPos(3)*(8/10);
    myPos(4) = myPos(4)*(8/10);
    set(gcf, 'Position', myPos);  % enlarge figure
    
    
    
    % TU Num timeline Panel
    subplot(xsize,ysize,posi.TimelinePanel)
    if ~isempty(plotTimeline)
        plot([0 0.5*(cnst.drawWhen:cnst.drawWhen:(cnst.drawWhen*numel(plotTimeline)))],...
            [4/3*pi*(r)^3 cell2mat(mysqueeze(plotTimeline.TU_Num))],'k','LineWidth',1.5);
    end
    axis([0 0.5*cnst.nSteps 0 700000])
    ylabel('tumor cells');
    title('tumor timeline');
    
    subplot(xsize,ysize,posi.minpHPanel);
    if ~isempty(plotTimeline)
        
        yyaxis left
        hold on
        
        plot([0 0.5*(cnst.drawWhen:cnst.drawWhen:(cnst.drawWhen*numel(plotTimeline)))],...
            [7.4 cell2mat(mysqueeze(plotTimeline.minpH))],'k','LineWidth',1.5);
        
        plot([0 0.5*(cnst.drawWhen:cnst.drawWhen:(cnst.drawWhen*numel(plotTimeline)))],...
            [5 cell2mat(mysqueeze(plotTimeline.minG))],'b','LineWidth',1.5);
        ylim([0 7.4])
        hold off
        ylabel('minimal pH/glucose')
        
        yyaxis right
        
        plot( [0 0.5*(cnst.drawWhen:cnst.drawWhen:(cnst.drawWhen*numel(plotTimeline)))],...
            [0.056 cell2mat(mysqueeze(plotTimeline.minO2))],'g','LineWidth',1.5);
        ylim([0 0.056])
    end
    xlim([0 0.5*cnst.nSteps ])
    ylabel('minimal oxygen')
    lgd=legend('pH', 'glucose', 'O2');
    lgd.Location='eastoutside';
    
    %
    %pH panel
    subplot(xsize, ysize, posi.spatialpHPanel)
    
    [X,Y,Z] = meshgrid(1:1:N1);
    data = -log10(mySystem.grid.ProtMap/1000);
    xslice = floor(N1/2);
    yslice = [];
    zslice = floor(N1/2);
    h=slice(X,Y,Z,data,xslice,yslice,zslice);
    set(h,'edgecolor','none');
    grid off
    colorbar;
    
    %PDL1 panel
    subplot(xsize, ysize, posi.PDLPanel)
    if ~isempty(plotTimeline)
        plot([0 0.5*(cnst.drawWhen:cnst.drawWhen:(cnst.drawWhen*numel(plotTimeline)))],...
            [0 cell2mat(mysqueeze(plotTimeline.PDL))],'k','LineWidth',1.5);
    end
    axis([0 0.5*cnst.nSteps 0 0.2])
    ylabel('Fraction of PDL1+ cells');
    
else
    if numel(mySystem.TU.TUcells)==0
        disp('no tumor cells, could not plot anything');
    end
end





end