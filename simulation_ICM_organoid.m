classdef simulation_ICM_organoid < handle
    % class simulation_ICM_organoid
    % |  A framework for the simulation of multicellular systems.
    % |
    % |  Example
    % |  -------
    % |  >> System = simulation_ICM_organoid()
    % |  >> for idx = time_steps
    % |  >>     System.TimeStep_Euler()
    % |  >> end
    % |  >> System.WriteData('final_time_point.mat')
    % |  
    % |  Notes
    % |  The class has been implemented to generate the results published
    % |  in Liebisch et al. (2020) BioRxiv (https://doi.org/10.1101/698928)
    
    
    properties (Access = public)
        % Discretisation
        dt      % time discretisation
        
        % Cell features
        N       % cell count
        X       % position data
        r       % cell size, given as radius
        r_max   % maximal cell size
        type    % cell fate (1: N-G-, 2: N+G+, 3: N+G-, 4: N-G+);
        
        % Mechanical forces
        F       % mechanical intercellular forces
        F1      % slope of the force potential
        F2      % multiplier for the force potential
        
        % Spatial parameter
        D       % Euklidean distances
        
        % Functions
        growth   % handle for the function defying cellular growth
        division % handle for the function defying division probabilities
        forcepot % handle for the applied forcepotential
        
        cell_cycle_time
        cell_cycle_time_list
        div_size_list
    end
    
    methods (Access = public)
        % User methods
        % |  User methods are called and accessible by the main programm. 
        % |  Using the constructor an instance of a simulation is generated 
        % |  (parameters are set within the constructor, please adapt if
        % |  required).
        % |  For testings a plot function is implemented.
        function self = simulation_ICM_organoid()
            % constructor simulation_ICM_organoid()
            % |  simulation_ICM_organoid()
            % |  
            % |  Declarations
            % |  ------------
            % |  N: integer
            % |      Initial cell count of the system
            % |  X: ndarray(N,3) of floats
            % |      Initial coordinates of the cells
            % |  r: ndarray(N,1) of floats
            % |      Initial cell size, given as their radii
            % |  r_max: float
            % |      Maximal cell size
            % |  F1: float
            % |      Slope of the forcepotential
            % |  F2: float
            % |      Multiplier for the forcepotential
            % |  growth: handle
            % |      Function defying cellular growth
            % |  division: handle
            % |      Function defying division probabilities
            % |  forcepot: handle
            % |      Function for the applied forcepotential
            self.dt = 0.001;
            self.N = 10;
            self.X = (rand(self.N,3)-0.5)*0.5;
            self.r = ones(self.N,1)*0.5;
            self.r_max = 1;
            self.F1 = 0.6;
            self.F2 = 1;
            
            self.growth = @(r, r_max) (r_max - r);
            
            self.division = @(r) 1 ./ (1 + exp(-(r - .8) / (1 / 16))) * 1.1 - 0.1;
            self.forcepot = @(d_is, d_sh, a) 2 * (exp(-2 * a * (d_is - d_sh)) - exp(-a * (d_is - d_sh)));
            
            self.type = ones(self.N,1);
        end
        function fTimeSetp_Euler(self)
            % function fTimeStep_Euler()
            % |  fTimeStep_Euler()
            % |
            % |  Notes
            % |  -----
            % |  Progresses through time for one timesetp.
            self.r = self.r + self.dt * self.growth(self.r, self.r_max);
            
            self.fDivision_H1() % _H2(); _H3(); _H4();
            self.fTriangulateDistanzes();
            self.fForcePotential();
            self.X = self.X + self.dt * self.F;
        end
        function fAssignCellFates(self, pNnGn, pNpGp, pNpGn, pNnGp)
            % function fAssignCellFates()
            % | fAssignCellFates()
            % | 
            % | Parameters
            % | ----------
            % | pNnGn: float
            % |     Probability of the cell Fata N-G- (0 <= pNnGn <= 1)
            % | pNpGp: float
            % |     Probability of the cell Fata N+G+ (0 <= pNpGp <= 1)
            % | pNpGn: float
            % |     Probability of the cell Fata N+G- (0 <= pNpGn <= 1)
            % | pNnGp: float
            % |     Probability of the cell Fata N-G+ (0 <= pNnGp <= 1)
            % | pNnGn + pNpGp + pNpGn + pNnGp = 1
            % | 
            % | Notes
            % | -----
            % | Assigns a cell fate to all cells in the system.
            for idx = 1:self.N
                rnd = rand();
                if rnd < pNnGn
                    self.type(idx) = 1;
                elseif rnd < pNpGp + pNnGn
                    self.type(idx) = 2;
                elseif rnd < pNpGn + pNpGp + pNnGn
                    self.type(idx) = 3;
                else % rnd < pNnGp + pNpGn + pNpGp + pNnGn
                    self.type(idx) = 4;
                end
            end
        end
    end
    methods (Access = public)
        % User methods
        % |  Here implemented functions are accessible by the main program 
        % |  in order to get current simulation variables (fGet[...]):
        % |  - N : cell count
        % |  - X : cell coordinates
        % |  - r : cell radii
        % |  - F : mechanical forces
        % |  The function fWriteData2csv(name) writes the variables X, r
        % |  and F to a .csv file at the location [name].
        % |  For testings a plot function is implemented. It is recommended
        % |  to use this function only for small simulations since it
        % |  increases the required computation time drastically.
        function N = fGetN(self)
            % function fGetN()
            % |  fGetN()
            % |  
            % |  Notes
            % |  -----
            % |  Returns the current cell count to the main program
            N = self.N;
        end
        function X = fGetX(self)
            % function fGetX()
            % |  fGetX()
            % |  
            % |  Notes
            % |  -----
            % |  Returns the current coordinates of the cells to the main 
            % |  program
            X = self.X;
        end
        function r = fGetR(self)
            % function fGetR()
            % |  fGetR()
            % |  
            % |  Notes
            % |  -----
            % |  Returns the radii of the cells to the main program
            r = self.r;
        end
        function t = fGetT(self)
            % function fGetT()
            % |  fGetT()
            % |  
            % |  Notes
            % |  -----
            % |  Returns the cell type of the cells to the main program
            t = self.type;
        end
        function F = fGetF(self)
            % function fGetF()
            % |  fGetF()
            % |  
            % |  Notes
            % |  -----
            % |  Returns the mechanical intercellular forces to the main
            % |  program
            F = self.F;
        end
        function D = fGetD(self)
            % function fGetD()
            % |  fGetD()
            % |  
            % |  Notes
            % |  -----
            % |  Returns the current distances between neighboured cells
            D = self.D;
        end
    end  
    methods (Access = private)
        % Privte methods
        % |  Here implemented functions are called by the class itself.
        % |  The functions are used to determine neighbouring cells and
        % |  to calculate their distances, to calculate the intercellular
        % |  mechanics and to determine which cells divide.
        function fTriangulateDistanzes(self)
            % function fTriangulateDistanzes()
            % | fTriangulateDistanzes()
            % | 
            % |  Return
            % |  ------
            % |  self.D: ndarray(self.N,self.N)
            % |      Euklidean distances between neighbouring cells. 
            % |      Distances between non neighbouring cells are set to 0.
            % |  self.Dx: ndarray(self.N,self.N)
            % |      Relative distances in x direction between neighbouring 
            % |      cells. Distances between non neighbouring cells are 
            % |      set to 0.
            % |  self.Dy: ndarray(self.N,self.N)
            % |      Relative distances in y direction between neighbouring 
            % |      cells. Distances between non neighbouring cells are 
            % |      set to 0.
            % |  self.Dz: ndarray(self.N,self.N)
            % |      Relative distances in z direction between neighbouring 
            % |      cells. Distances between non neighbouring cells are 
            % |      set to 0.
            self.D  = zeros(self.N);
            
            DTR = delaunayTriangulation(self.X(:,1),self.X(:,2),self.X(:,3));
            dtr = DTR.ConnectivityList;
            dim = size(dtr);
                
            for idx1 = 1:dim(2)
                for idx2 = idx1:dim(2)
                    if idx1 ~= idx2
                        n1 = dtr(:,idx1);
                        n2 = dtr(:,idx2);
                        
                        lindx1=n1+self.N*(n2-1);
                        lindx2=n2+self.N*(n1-1);
                        
                        dd = ((self.X(n1,1)-self.X(n2,1)).^2 ...
                             +(self.X(n1,2)-self.X(n2,2)).^2 ...
                             +(self.X(n1,3)-self.X(n2,3)).^2).^0.5;
                        self.D(lindx1)  = dd;
                        self.D(lindx2)  = dd;                    
                    end
                end
            end  
        end
        function fForcePotential(self)
            % function fForcePotential()
            % | fForcePotential()
            % | 
            % |  Return
            % |  ------
            % |  self.F: ndarray(self.N,3)
            % |      Mechanical intercellular forces.

            self.F = zeros(self.N,3);
            for idx = 1:self.N
                neighbours = find(self.D(idx,:)>0);
                if ~isempty(neighbours)
                    D_is = self.D(idx,neighbours)';
                    D_sh = self.r(idx)+self.r(neighbours);
                    
                    REP = self.forcepot(D_is,D_sh,self.F1);
                    
                    self.F(idx,:) = self.F2 * (REP' * (self.X(idx,:) - self.X(neighbours,:)));
                end
            end
        end
        function fDivision_H1(self)
            % function fDivision()
            % | fDivision()
            % | 
            % |  Return
            % |  ------
            % |  self.r (new): ndarray(N,1) of floats
            % |      Radii of the cells
            % |  self.N (new): integer
            % |      Cell count
            % |  self.X (new): ndarray(N,3) of floats
            % |      Coordinates of the cells   
            % | 
            % |  Notes
            % |  -----
            % |   Resmbles hypothesis 1 of cell fate heredity i.e. strict
            % |   cell fate heredity, no cell fate switches alllowed.   
            p_division = self.dt * self.division(self.r);
            div = find(p_division > rand(self.N, 1));
            if ~isempty(div)
                l = length(div);
                
                self.N = self.N + l;
                
                division_plane = rand(l,3)-0.5;
                division_plane = division_plane ./ sum(division_plane.^2, 2) .^ 0.5;
                
                old_r = self.r(div);
                new_r = 0.7937 * self.r(div); % (1/2)^(1/3) = 0.7937
                self.r(div) = new_r ;
                self.r = cat(1, self.r, new_r);
                
                self.type = cat(1, self.type, self.type(div));
                
                old_X = self.X(div, :);
                self.X(div, :) = old_X - (old_r - new_r) .* division_plane;
                X_new = old_X + (old_r - new_r) .* division_plane;
                self.X = cat(1, self.X, X_new);
            end
        end
        function fDivision_H2(self)
            % function fDivision()
            % | fDivision()
            % | 
            % |  Return
            % |  ------
            % |  self.r (new): ndarray(N,1) of floats
            % |      Radii of the cells
            % |  self.N (new): integer
            % |      Cell count
            % |  self.X (new): ndarray(N,3) of floats
            % |      Coordinates of the cells   
            % | 
            % |  Notes
            % |  -----
            % |   Resmbles hypothesis 2 of cell fate heredity i.e. cell
            % |   fate switches from N+G+ to N-G+ are allowed.    
            p_division = self.dt * self.division(self.r);
            div = find(p_division > rand(self.N, 1));
            if ~isempty(div)
                l = length(div);
                
                self.N = self.N + l;
                
                division_plane = rand(l, 3) - 0.5;
                division_plane = division_plane ./ sum(division_plane.^2, 2) .^ 0.5;
                
                old_r = self.r(div);
                new_r = 0.7937 * self.r(div); % (1/2)^(1/3) = 0.7937
                self.r(div) = new_r ;
                self.r = cat(1, self.r, new_r);
                
                new_type = self.type(div);
                
                for idx = 1:l
                    if self.type(div(idx)) == 2
                        rnd = rand();
                        if rnd < 0.5
                            self.type(div(idx)) = 4;
                        end
                        rnd = rand();
                        if rnd < 0.5
                            new_type(idx) = 4;
                        end
                    end
                end
                
                self.type = cat(1, self.type, new_type);
                
                old_X = self.X(div, :);
                self.X(div, :) = old_X - (old_r - new_r) .* division_plane;
                X_new = old_X + (old_r - new_r) .* division_plane;
                self.X = cat(1, self.X, X_new);
            end
        end
        function fDivision_H3(self)
            % function fDivision()
            % | fDivision()
            % | 
            % |  Return
            % |  ------
            % |  self.r (new): ndarray(N,1) of floats
            % |      Radii of the cells
            % |  self.N (new): integer
            % |      Cell count
            % |  self.X (new): ndarray(N,3) of floats
            % |      Coordinates of the cells   
            % | 
            % |  Notes
            % |  -----
            % |   Resmbles hypothesis 3 of cell fate heredity i.e. cell
            % |   fate switches from N-G- to N+G+, N+G+ to N+G- or N-G+, 
            % |   N+G- to N-G+ are allowed.    
            p_division = self.dt * self.division(self.r);
            div = find(p_division > rand(self.N, 1));
            if ~isempty(div)
                l = length(div);
                
                self.N = self.N + l;
                
                division_plane = rand(l, 3) - 0.5;
                division_plane = division_plane ./ sum(division_plane.^2, 2) .^ 0.5;
                
                old_r = self.r(div);
                new_r = 0.7937 * self.r(div); % (1/2)^(1/3) = 0.7937
                self.r(div) = new_r ;
                self.r = cat(1, self.r, new_r);
                
                new_type = self.type(div);
                
                for idx = 1:l
                    if self.type(div(idx)) == 1
                        rnd = rand();
                        if rnd < 0.1061
                            self.type(div(idx)) = 2;
                        end
                        rnd = rand();
                        if rnd < 0.1061
                            new_type(idx) = 2;
                        end
                    elseif self.type(div(idx)) == 2
                        rnd = rand();
                        if rnd < 0.3932
                            self.type(div(idx)) = 3;
                        elseif rnd < 0.3932 * 2
                            self.type(div(idx)) = 4;
                        end
                        rnd = rand();
                        if rnd < 0.3932
                            new_type(idx) = 3;
                        elseif rnd < 0.3932 * 2
                            new_type(idx) = 4;
                        end
                    elseif self.type(div(idx)) == 3
                        rnd = rand();
                        if rnd < 0.1498
                            self.type(div(idx)) = 4; 
                        end
                        rnd = rand();
                        if rnd < 0.1498
                            new_type(idx) = 4; 
                        end
                    end
                end
                
                self.type = cat(1, self.type, new_type);
                
                old_X = self.X(div, :);
                self.X(div, :) = old_X - (old_r - new_r) .* division_plane;
                X_new = old_X + (old_r - new_r) .* division_plane;
                self.X = cat(1, self.X, X_new);
            end
        end
        function fDivision_H4(self)
            % function fDivision()
            % | fDivision()
            % | 
            % |  Return
            % |  ------
            % |  self.r (new): ndarray(N,1) of floats
            % |      Radii of the cells
            % |  self.N (new): integer
            % |      Cell count
            % |  self.X (new): ndarray(N,3) of floats
            % |      Coordinates of the cells   
            % | 
            % |  Notes
            % |  -----
            % |   Resmbles hypothesis 4 of cell fate heredity i.e. cell
            % |   fate switches from N-G- to N+G+, N+G+ to N+G- or N-G+, 
            % |   N+G- to N-G+, N-G+ to N+G- are allowed.    
            p_division = self.dt * self.division(self.r);
            div = find(p_division > rand(self.N, 1));
            if ~isempty(div)
                l = length(div);
                
                self.N = self.N + l;
                
                division_plane = rand(l, 3) - 0.5;
                division_plane = division_plane ./ sum(division_plane.^2, 2) .^ 0.5;
                
                old_r = self.r(div);
                new_r = 0.7937 * self.r(div); % (1/2)^(1/3) = 0.7937
                self.r(div) = new_r ;
                self.r = cat(1, self.r, new_r);
                
                new_type = self.type(div);
                
                for idx = 1:l
                    if self.type(div(idx)) == 1
                        rnd = rand();
                        if rnd < 0.1061
                            self.type(div(idx)) = 2;
                        end
                        rnd = rand();
                        if rnd < 0.1061
                            new_type(idx) = 2;
                        end
                    elseif self.type(div(idx)) == 2
                        rnd = rand();
                        if rnd < 0.3932
                            self.type(div(idx)) = 3;
                        elseif rnd < 0.3932 * 2
                            self.type(div(idx)) = 4;
                        end
                        rnd = rand();
                        if rnd < 0.3932
                            new_type(idx) = 3;
                        elseif rnd < 0.3932 * 2
                            new_type(idx) = 4;
                        end
                    elseif self.type(div(idx)) == 3
                        rnd = rand();
                        if rnd < 0.2196
                            self.type(div(idx)) = 4; 
                        end
                        rnd = rand();
                        if rnd < 0.2196
                            new_type(idx) = 4; 
                        end
                    elseif self.type(div(idx)) == 4
                        rnd = rand();
                        if rnd < 0.0750
                            self.type(div(idx)) = 3; 
                        end
                        rnd = rand();
                        if rnd < 0.0750
                            new_type(idx) = 3; 
                        end
                    end
                end
                
                self.type = cat(1, self.type, new_type);
                
                old_X = self.X(div, :);
                self.X(div, :) = old_X - (old_r - new_r) .* division_plane;
                X_new = old_X + (old_r - new_r) .* division_plane;
                self.X = cat(1, self.X, X_new);
            end
        end
    end
end