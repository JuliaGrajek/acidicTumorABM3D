%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef SIMengine_interface < handle
    properties (SetAccess = public, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = SIMengine_interface(varargin)
            this.objectHandle = SIMengine_interface_mex('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            SIMengine_interface_mex('delete', this.objectHandle);
        end
        
        function varargout = initialize(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('initialize', this.objectHandle, varargin{:});
        end
        
        function varargout = initializeFromState(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('initializeFromState', this.objectHandle, varargin{:});
        end
        
        
        function varargout = getState(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('getState', this.objectHandle, varargin{:});
        end
        
        function varargout = readState(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('readState', this.objectHandle, varargin{:});
        end
        
        function TU_go_grow_die(this, varargin)
            SIMengine_interface_mex('TU_go_grow_die', this.objectHandle, varargin{:});
        end
        
        function TU_just_die(this, varargin)
            SIMengine_interface_mex('TU_just_die', this.objectHandle, varargin{:});
        end
        
        function modulateIFNgMap(this, varargin)
            SIMengine_interface_mex('modulateIFNgMap', this.objectHandle, varargin{:});
        end
        
        function decayIFNgMap(this, varargin)
            SIMengine_interface_mex('decayIFNgMap', this.objectHandle, varargin{:});
        end
        
        function modulateHypoxia(this, varargin)
            SIMengine_interface_mex('modulateHypoxia', this.objectHandle, varargin{:});
        end
        
        function IMinflux(this, varargin)
            SIMengine_interface_mex('IMinflux', this.objectHandle, varargin{:});
        end

        
        function lymphocytesAct(this, varargin)
            SIMengine_interface_mex('lymphocytesAct', this.objectHandle, varargin{:});
        end
        
        
        function updateNecroMap(this, varargin)
            SIMengine_interface_mex('updateNecroMap', this.objectHandle, varargin{:});
        end
        
        function updateGlucMap(this, varargin)
            SIMengine_interface_mex('updateGlucMap', this.objectHandle, varargin{:});
        end
        
        function updateProtMap(this, varargin)
            SIMengine_interface_mex('updateProtMap', this.objectHandle, varargin{:});
        end
        
        function updateChemoMap(this, varargin)
            SIMengine_interface_mex('updateChemoMap', this.objectHandle, varargin{:});
        end
        
        function seedFibrosis(this, varargin)
            SIMengine_interface_mex('seedFibrosis', this.objectHandle, varargin{:});
        end
        
        function clearProtonSources(this, varargin)
            SIMengine_interface_mex('clearProtonSources', this.objectHandle, varargin{:});
        end
		
		function varargout = TUcellsNum(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('TUcellsNum', this.objectHandle, varargin{:});
        end
        
        
    end
end