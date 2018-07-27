function [y, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x, RMSin)
% Version date: 2017-3-10
%function [y, RMSout, Idc, Vdc] = IMS2015_PA_meas_Pin(x, RMSin)
%Executes measurement of remote PA using IQ-data x at power level RMSin. x
%is sampled at 200 MHz. Please note: the bandwidth is limited to 160 MHz.
%This means you should not put any useful signal components outside of the
%range [-80, 80] MHz because it will be distorted in the signal generation
%process. You can also not use the returned data outside of the range
%[-80, 80] MHz because it is filtered by the signal analyzer.
%
%Returns IQ-data y with RMS power level RMSout (into a 50 Ohm load). y is
%scaled such that y it is given in Volts. y is sampled at 200 MHz and has a
%bandwidth of 160 MHz.
%
%There are limitations on the rms and peak power levels output from the
%signal generator: maximum rms output depends on the peak-to-average ratio 
%and the maximum peak output is -11 dBm. These limits are put there to ensure 
%that the PA is not unintentionally damaged by excessive power levels.
%
%Inputs:
% Note that the client can be run in two modes. Firstly by specifying both x 
% and RMSin. Or, only by specificying x,which is given in Volts (into a 50 Ohm
%    load). The firstcase is shortly described below, after which an example is given.
%
% x - input IQ-data sampled at 200 MHz. Note that x is always rescaled to
%    have max(abs(x))==1 to utilize the DAC as well as possible. x has a
%    limited peak-to-average power ratio (PAPR) of 20 dB.
%    Maximum 1 000 000 samples are allowed. Longer files are cut to 1 000 000
%    samples at the server.
%    Note that x must have at least 1000 samples. 
%    Note that x should have an even number of samples.
% RMSin - rms power level setting in dBm of the signal generator. The
%        approximate gain of the amplifier stages preceding the PA is 40
%        dB. RMSin is limited to( PEAKin - PAPR(x) ) dBm rms and the 
%        maximum peak power is limited to -11 dBm.
%
% Usage examples:
%   - for x as a voltage mode ->  [y, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1(x);  
%   - for RMSin mode          ->  [y, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1(x,RMSin); 
%
%Outputs:
%y - measured IQ-data sampled at 200 MHz with a bandwidth of 160 MHz. y is
%    always given in Volts. 
%    y has the same number of samples as x.
%    The sampling is started in a completely random fashion so you need to
%    do a synchronization before using the output data.
%RMSout - rms output power of y in dBm.
%Idc - measured current
%Vdc - voltage (20V)
%

warning off backtrace
if nargin < 2
   RMSin = 10*log10( norm(x)^2/50/length(x)) + 30;                         %convert V-rms to dBm
end

[status,x] = check_inputs_WebLab(x, RMSin);                                %Input validation check
                                                                           
if status  
    GA_call();                                                             %Google Analytics
                                                                           %If everything is ok, proceed
    [InData]= createInCluster_WebLab(x,RMSin);                             %Create system specific Cluster struc
    [FormedCluster] = form_data_cluster(InData);                           %Convert input cluster to bin and apply hashing
    [h,h_ok] = filewrite(FormedCluster);                                   %write data to server, get handle back

    maxTry = 100;                                                          %try for maxTry*pauselen seconds to download data, otherwise return empty y
    pauselen = 3;
    
    [OutData]= fileread(h,h_ok, maxTry, pauselen);
    
    if OutData.status 
    [y,RMSout,Idc,Vdc]= readOutData_WebLab(OutData);
    disp('Results are ready')
    else   
    y = []; RMSout = []; Idc = []; Vdc = [];
    end
    
else
    error('Power check failed. Rerun or reduce the RMS Power Level')
    y = []; RMSout = []; Idc = []; Vdc = [];
end

end %end of main function

%% Sub Functions for specific system

function [status,x] = check_inputs_WebLab(x, RMSin)

status=true;

if isempty(x)
    warning('Invalid IQ Data. Empty Data Array')
    status = false;return     
elseif ~all(isfinite(x))
    warning('Invalid IQ Data. Data contain NaN or Inf value(s)')
    status = false;return 
elseif ~isvector(x)
    warning('Invalid IQ Data. IQ data must be a vector')
    status = false;return 
elseif length(x)>1e6
    warning('x is too long, cutting to 1 000 000 samples.')
    x = x(1:1e6);
elseif length(x)<1000
    warning('x must be more than 1000 samples.')
    status = false;return 
elseif isreal(x)
    x= complex(x,0);
end

if isempty(RMSin)
    warning('Invalid RMSin Data. Empty Data value')
    status = false;return     
elseif ~all(size(RMSin)==[1 1])
    warning('Invalid RMSin value. Single value is required')
    status = false;return 
elseif ~isfinite(RMSin)
    warning('Invalid RMSin value. Data contain NaN or Inf value(s)')
    status = false;return     
end


x = x(:); 
% Normalization of data
PeakReal=max(abs(real(x)));  
PeakImag=max(abs(imag(x)));
Peak=max([PeakReal PeakImag]);
if Peak ~=0
    x=double(x./Peak);
else
    x=double(x);
end

PAPR = 20*log10(max(abs(x))*sqrt(length(x))/norm(x)); 

if RMSin > -8 - PAPR*0.77 
    warning('rms input power too high.')
    status = false;return  
elseif PAPR > 20
    warning('Signal PAPR too high.')
    status = false;return 
elseif RMSin+PAPR > -8
    warning('Peak input power too high.')
    status = false;return 
end

end

function [InData]= createInCluster_WebLab(x,RMSin)
%Create the in cluster with the settings and data 
%Data struct creation
InData.Client_version=single(1.1);
InData.OpMode        ='A5Z6UNud';
InData.PowerLevel    =RMSin;
InData.SA_Samples    =length(x);
InData.Re_data       =real(x);
InData.Im_data       =imag(x);

end

function [y,RMSout,Idc,Vdc]= readOutData_WebLab(OutData)

y= complex(OutData.b3_re,OutData.b3_im);
RMSout=OutData.RMS_out;
Idc=OutData.DCmeas.Id;
Vdc=OutData.DCmeas.Vd;

end

%% Sub Functions

function [OutData] = fileread(h,h_ok, maxTry, pauselen)

%Read data from the server

ind1 = strfind(h, 'out');
hk= [h(ind1:end-3) 'dat'];

for k=1:maxTry
    pause(pauselen) 
    try
        
        tmp_nam_ok = tempname; %temporary filename in case reading success
        urlwrite(h_ok,tmp_nam_ok);
        delete(tmp_nam_ok)
        
        tmp_nam = tempname; %temporary filename in case reading success
        urlwrite(h,tmp_nam);
        unzip(tmp_nam);
        
        fid = fopen(hk,'r','l');
        OutCluster_bin = fread(fid,'*uint8'); %read outdata binary
        fclose(fid); %close
        delete(hk) %delete temporary data
        
        OutData = deserialize(OutCluster_bin);   %convert outdata binary to stuct
        
        if isfield(OutData, 'message')
            disp(OutData.message);
        end
       
        error_handle(OutData);
        timeout=false;
        break;
       
    catch
        OutData.status=false;
        timeout=true;
    end
   
end

if timeout
   warning('Connection timeout. Check your internet connection and try again.') 
end

end %end of fileread

function error_handle(OutData)

if OutData.error_hanlde==int32(-999000);
    %No warning 
elseif OutData.error_hanlde==int32(-999001);
    warning('Measurement System failure. Please repeat the measurement.')    
elseif OutData.error_hanlde==int32(-999002);
    warning('Invalid input data. Please verify your settings.')    
elseif OutData.error_hanlde==int32(-999003);
    warning('You are using an old Client version. Please update.')
elseif OutData.error_hanlde==int32(-999004);
    warning('Corrupted data transfer. Please repeat the measurement.') 
end

end

function [s,s_ok, hname] = filewrite(FormedCluster)
%%Write formatted InputCluster to server, get a handle to the measured data back

%Use URLREADPOST from Dan Ellis to upload data, please note license
%agreement of URLREADPOST
urlChar = 'http://dpdcompetition.com/rfweblab/matlab/upload.php';
params = {'myFile',FormedCluster};
out = urlreadpost(urlChar,params);

%parse output file name
ind1 = strfind(out, 'output_'); 
ind1 = ind1(1);
%parse output ok file name
ind2 = strfind(out,'.dat'); %file always ends with .dat, otherwise not valid file
ind2 = ind2(1);

fname = out(ind1:ind2+3); %this is the file name on the server
fnameok = ['ok_' out(ind1+7:ind2+3)];
urlChar2 = 'http://dpdcompetition.com/rfweblab/matlab/files/'; 
s = [urlChar2 fname(1:end-3) 'zip']; %URL to where the data is located
s_ok= [urlChar2 fnameok];
hname = fname(7:end); %only last part of file name

disp('Waiting for measurement...')
end %end of filewrite

function [output,status] = urlreadpost(urlChar,params)
%URLREADPOST Returns the contents of a URL POST method as a string.
%   S = URLREADPOST('URL',PARAMS) passes information to the server as
%   a POST request.  PARAMS is a cell array of param/value pairs.
%   
%   Unlike stock urlread, this version uses the multipart/form-data
%   encoding, and can thus post file content.  File data is 
%   encoded as a value element of numerical type (e.g. uint8)
%   in PARAMS.  For example:
%
%   f = fopen('music.mp3');
%   d = fread(f,Inf,'*uint8');  % Read in byte stream of MP3 file
%   fclose(f);
%   str = urlreadpost('http://developer.echonest.com/api/upload', ...
%           {'file',d,'version','3','api_key','API-KEY','wait','Y'});
%
%   ... will upload the mp3 file to the Echo Nest Analyze service.
%
%  Based on TMW's URLREAD.  Note that unlike URLREAD, there is no
%  METHOD argument
%  2010-04-07 Dan Ellis dpwe@ee.columbia.edu
% This function requires Java.
if ~usejava('jvm')
   error('MATLAB:urlreadpost:NoJvm','URLREADPOST requires Java.');
end
import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;
% Be sure the proxy settings are set.
com.mathworks.mlwidgets.html.HTMLPrefs.setProxySettings
% Check number of inputs and outputs.
narginchk(2,2)
nargoutchk(0,2)
if ~ischar(urlChar)
    error('MATLAB:urlreadpost:InvalidInput','The first input, the URL, must be a character array.');
end
% Do we want to throw errors or catch them?
if nargout == 2
    catchErrors = true;
else
    catchErrors = false;
end
% Set default outputs.
output = '';
status = 0;
% Create a urlConnection.
[urlConnection,errorid,errormsg] = urlreadwrite_old(mfilename,urlChar);
if isempty(urlConnection)
    if catchErrors, return
    else error(errorid,errormsg);
    end
end
% POST method.  Write param/values to server.
% Modified for multipart/form-data 2010-04-06 dpwe@ee.columbia.edu
%    try
        urlConnection.setDoOutput(true);
        boundary = '***********************';
        urlConnection.setRequestProperty( ...
            'Content-Type',['multipart/form-data; boundary=',boundary]);
        printStream = java.io.PrintStream(urlConnection.getOutputStream);
        % also create a binary stream
        dataOutputStream = java.io.DataOutputStream(urlConnection.getOutputStream);
        eol = [char(13),char(10)];
        for i=1:2:length(params)
          printStream.print(['--',boundary,eol]);
          printStream.print(['Content-Disposition: form-data; name="',params{i},'"']);
          if ~ischar(params{i+1}) 
            % binary data is uploaded as an octet stream
            % Echo Nest API demands a filename in this case
            printStream.print(['; filename="dummy.dat"',eol]);
            printStream.print(['Content-Type: application/octet-stream',eol]);
            printStream.print([eol]);
            dataOutputStream.write(params{i+1},0,length(params{i+1}));
            printStream.print([eol]);
          else
            printStream.print([eol]);
            printStream.print([eol]);
            printStream.print([params{i+1},eol]);
          end
        end
        printStream.print(['--',boundary,'--',eol]);
        printStream.close;
%    catch
%        if catchErrors, return
%        else error('MATLAB:urlread:ConnectionFailed','Could not POST to URL.');
%        end
%    end

% Read the data from the connection.
try
    inputStream = urlConnection.getInputStream;
    byteArrayOutputStream = java.io.ByteArrayOutputStream;
    % This StreamCopier is unsupported and may change at any time.
    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
    isc.copyStream(inputStream,byteArrayOutputStream);
    inputStream.close;
    byteArrayOutputStream.close;
    output = native2unicode(typecast(byteArrayOutputStream.toByteArray','uint8'),'UTF-8');
catch
    if catchErrors, return
    else error('MATLAB:urlreadpost:ConnectionFailed','Error downloading URL. Your network connection may be down or your proxy settings improperly configured.');
    end
end
status = 1;
end %end of URLREADPOST

function [urlConnection,errorid,errormsg] = urlreadwrite_old(fcn,urlChar)
%URLREADWRITE A helper function for URLREAD and URLWRITE.

%   Matthew J. Simoneau, June 2005
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3.6.1 $ $Date: 2009/01/30 22:37:42 $

% Default output arguments.
urlConnection = [];
errorid = '';
errormsg = '';

% Determine the protocol (before the ":").
protocol = urlChar(1:min(find(urlChar==':'))-1);

% Try to use the native handler, not the ice.* classes.
switch protocol
    case 'http'
        try
            handler = sun.net.www.protocol.http.Handler;
        catch exception %#ok
            handler = [];
        end
    case 'https'
        try
            handler = sun.net.www.protocol.https.Handler;
        catch exception %#ok
            handler = [];
        end
    otherwise
        handler = [];
end

% Create the URL object.
try
    if isempty(handler)
        url = java.net.URL(urlChar);
    else
        url = java.net.URL([],urlChar,handler);
    end
catch exception %#ok
    errorid = ['MATLAB:' fcn ':InvalidUrl'];
    errormsg = 'Either this URL could not be parsed or the protocol is not supported.';
    return
end

% Get the proxy information using MathWorks facilities for unified proxy
% prefence settings.
mwtcp = com.mathworks.net.transport.MWTransportClientPropertiesFactory.create();
proxy = mwtcp.getProxy(); 


% Open a connection to the URL.
if isempty(proxy)
    urlConnection = url.openConnection;
else
    urlConnection = url.openConnection(proxy);
end
end %end of urlreadwrite_old

function GA_call()
URL = 'http://www.google-analytics.com/__utm.gif';
utmwv='5.4.7'; %tracking version
utms='3'; %visit number, must probably match with some cookie setting later
utmn='1378994324'; %site hash
utmhn='dpdcompetition.com'; %site
utmcs='UTF-8'; %encoding
utmsr='1680x1050'; %screen res
utmvp='1664x486'; %browser res
utmsc='24-bit'; %color depth
utmul='en-us'; %language
utmje='1'; %JS enabled
utmfl='12.0 r0'; %Flash version
utmdt='Power Amplifier Linearization through Digital Predistortion (DPD) | IMS2014 Student Design Competition'; %title of page
utmhid='1518660351'; %random number for AdSense
utmr='0'; %referral, complete url
utmp='/'; %page request of the current page
utmht='1391522837321'; %some kind of time stamp
utmac='UA-47651321-1'; %account string
%following is a cookie report
utmcc='__utma=26860264.1502292211.1391522828.1391522828.1391522828.1;+__utmz=26860264.1391522828.1.1.utmcsr=(direct)|utmccn=(direct)|utmcmd=(none);';
utmu='q~'; %some info setting for google
params = {'utmwv',utmwv, 'utms',utms, 'utmn',utmn, 'utmhn',utmhn, 'utmcs',utmcs,...
    'utmsr',utmsr, 'utmvp',utmvp, 'utmsc',utmsc, 'utmul',utmul, 'utmje',utmje,...
    'utmfl',utmfl, 'utmdt',utmdt, 'utmhid',utmhid, 'utmr',utmr, 'utmp',utmp, 'utmht',utmht,...
    'utmac',utmac, 'utmcc',utmcc, 'utmu',utmu};
s=urlread(URL,'get',params); %make call to gif with all the parameters using GET method
 
end %end of GA_call

%% Aux Functions
function [FormedCluster] = form_data_cluster(InData)
InCluster_bin = serialize(InData);
Opt = struct('Format', 'uint8', 'Method', 'MD5');
hash=dataHash(InCluster_bin, Opt).';
FormedCluster=[hash;InCluster_bin];




function m = serialize(v)
% SERIALIZE converts a matlab object into a compact (but uncompressed)
% series of bytes.
%
% m = SERIALIZE(v)
%
% v is a matlab object. It can be any combination of
% structs, cells, and arrays. Other object types are not supported.
% Note that numeric types in matlab are double by default. To save space
% you might want to convert them to integers like this:
%
% foo = int16([1 2 3]);
%
% They are automatically converted back to doubles by DESERIALIZE.
%
% Limitations: No object can have more than 255 dimensions, and each
% dimension must be smaller than 2^32. Also structs cannot have more
% than 255 fields.
%
% EXAMPLE:
%
% v(1).red = {'blue', uint16([1 2 3])};
% v(2).green = [4.5 6.7 8.9];
% m = serialize(v);
% v2 = deserialize(m);
%
% By Tim Hutt, 19/11/2010
%
% Updated 16/12/2010 - Fix bug with matrices.

if isnumeric(v) || islogical(v) || ischar(v) % Matrix type thing.
    m = serializeMatrix(v);
elseif isstruct(v)
    m = serializeStruct(v);
elseif iscell(v)
    m = serializeCell(v);
else
    error('Unknown class');
end

end

function m = serializeMatrix(v)
m = uint8([]);
% Data type.
m = [m; classToByte(class(v))];

% Number of dimensions.
m = [m; ndims(v)];

% Dimensions.
for ii = 1:ndims(v)
    m = [m; typecast(uint32(size(v, ii)), 'uint8').'];
end

% Data.
if ischar(v)
    m = [m; uint8(v(:))];
else
    m = [m; typecast(v(:).', 'uint8').'];
end
end

function b = classToByte(cls)
switch cls
    case 'double'
        b = 0;
    case 'single'
        b = 1;
    case 'logical'
        b = 2;
    case 'char'
        b = 3;
    case 'int8'
        b = 4;
    case 'uint8'
        b = 5;
    case 'int16'
        b = 6;
    case 'uint16'
        b = 7;
    case 'int32'
        b = 8;
    case 'uint32'
        b = 9;
    case 'int64'
        b = 10;
    case 'uint64'
        b = 11;
    otherwise
        error('Unknown class');
end
end


function m = serializeCell(v)
m = uint8([]);
% Data type.
m = [m; 254]; % 254 = cell.

% Number of dimensions.
m = [m; ndims(v)];

% Dimensions.
for ii = 1:ndims(v)
    m = [m; typecast(uint32(size(v, ii)), 'uint8').'];
end

% Just serialize each member.
for ii = 1:numel(v)
    m = [m; serialize(v{ii})];
end
end

% Struct array. A plain struct is just a struct array of size 1.
function m = serializeStruct(v)
m = uint8([]);
% Data type.
m = [m; 255]; % 255 = struct.


% Field names.
fieldNames = fieldnames(v);

if numel(fieldNames) > 255
    error('Too many fields!');
end

% Number of field names.
m = [m; numel(fieldNames)];

for ii = 1:numel(fieldNames)
    m = [m; typecast(uint32(numel(fieldNames{ii})), 'uint8').'; uint8(fieldNames{ii}(:))];
end


% Number of dimensions.
m = [m; ndims(v)];

% Dimensions.
for ii = 1:ndims(v)
    m = [m; typecast(uint32(size(v, ii)), 'uint8').'];
end

% Now for go through each one.
% This is slightly redundant because we encode the type info lots of
% times, but it is only one byte so meh.
for ii = 1:numel(v)
    for ff = 1:numel(fieldNames)
        m = [m; serialize(v(ii).(fieldNames{ff}))];
    end
end
end



function Hash = dataHash(Data, Opt)
% DATAHASH - Checksum for Matlab array of any type
% This function creates a hash value for an input of any type. The type and
% dimensions of the input are considered as default, such that UINT8([0,0]) and
% UINT16(0) have different hash values. Nested STRUCTs and CELLs are parsed
% recursively.
%
% Hash = DataHash(Data, Opt)
% INPUT:
%   Data: Array of these built-in types:
%           (U)INT8/16/32/64, SINGLE, DOUBLE, (real/complex, full/sparse)
%           CHAR, LOGICAL, CELL (nested), STRUCT (scalar or array, nested),
%           function_handle.
%   Opt:  Struct to specify the hashing algorithm and the output format.
%         Opt and all its fields are optional.
%         Opt.Method: String, known methods for Java 1.6 (Matlab 2011b):
%              'SHA-1', 'SHA-256', 'SHA-384', 'SHA-512', 'MD2', 'MD5'.
%            Call DataHash without inputs to get a list of available methods.
%            Default: 'MD5'.
%         Opt.Format: String specifying the output format:
%            'hex', 'HEX':      Lower/uppercase hexadecimal string.
%            'double', 'uint8': Numerical vector.
%            'base64':          Base64 encoded string, only printable ASCII
%                               characters, shorter than 'hex', no padding.
%            Default: 'hex'.
%         Opt.Input: Type of the input as string, not case-sensitive:
%            'array': The contents, type and size of the input [Data] are
%                     considered  for the creation of the hash. Nested CELLs
%                     and STRUCT arrays are parsed recursively. Empty arrays of
%                     different type reply different hashs.
%            'file':  [Data] is treated as file name and the hash is calculated
%                     for the files contents.
%            'bin':   [Data] is a numerical, LOGICAL or CHAR array. Only the
%                     binary contents of the array is considered, such that
%                     e.g. empty arrays of different type reply the same hash.
%            'ascii': Same as 'bin', but only the 8-bit ASCII part of the 16-bit
%                     Matlab CHARs is considered.
%            Default: 'array'.
%
% OUTPUT:
%   Hash: String, DOUBLE or UINT8 vector. The length depends on the hashing
%         method.
%
% EXAMPLES:
% % Default: MD5, hex:
%   DataHash([])                % 5b302b7b2099a97ba2a276640a192485
% % MD5, Base64:
%   Opt = struct('Format', 'base64', 'Method', 'MD5');
%   DataHash(int32(1:10), Opt)  % +tJN9yeF89h3jOFNN55XLg
% % SHA-1, Base64:
%   S.a = uint8([]);
%   S.b = {{1:10}, struct('q', uint64(415))};
%   Opt.Method = 'SHA-1';
%   Opt.Format = 'HEX';
%   DataHash(S, Opt)            % 18672BE876463B25214CA9241B3C79CC926F3093
% % SHA-1 of binary values:
%   Opt = struct('Method', 'SHA-1', 'Input', 'bin');
%   DataHash(1:8, Opt)          % 826cf9d3a5d74bbe415e97d4cecf03f445f69225
% % SHA-256, consider ASCII part only (Matlab's CHAR has 16 bits!):
%   Opt.Method = 'SHA-256';
%   Opt.Input  = 'ascii';
%   DataHash('abc', Opt)
%       % ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad
%   % Or equivalently:
%   Opt.Input = 'bin';
%   DataHash(uint8('abc'), Opt)
%
% NOTES:
%   Function handles and user-defined objects cannot be converted uniquely:
%   - The subfunction ConvertFuncHandle uses the built-in function FUNCTIONS,
%     but the replied struct can depend on the Matlab version.
%   - It is tried to convert objects to UINT8 streams in the subfunction
%     ConvertObject. A conversion by STRUCT() might be more appropriate.
%   Adjust these subfunctions on demand.
%
%   MATLAB CHARs have 16 bits! Use Opt.Input='ascii' for comparisons with e.g.
%   online hash generators.
%
%   Matt Raum suggested this for e.g. user-defined objects:
%     DataHash(getByteStreamFromArray(Data)
%   This works very well, but unfortunately getByteStreamFromArray is
%   undocumented, such that it might vanish in the future or reply different
%   output.
%
%   For arrays the calculated hash value might be changed in new versions.
%   Calling this function without inputs replies the version of the hash.
%
%   The C-Mex function GetMD5 is 2 to 100 times faster, but obtains MD5 only:
%   http://www.mathworks.com/matlabcentral/fileexchange/25921
%
% Tested: Matlab 7.7, 7.8, 7.13, 8.6, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2011-2016 matlab.2010(a)n(MINUS)simon.de
%
% See also: TYPECAST, CAST.
%
% Michael Kleder, "Compute Hash", no structs and cells:
%   http://www.mathworks.com/matlabcentral/fileexchange/8944
% Tim, "Serialize/Deserialize", converts structs and cells to a byte stream:
%   http://www.mathworks.com/matlabcentral/fileexchange/29457

% $JRev: R-H V:033 Sum:R+m7rAPNLvlw Date:18-Jun-2016 14:33:17 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\GLFile\DataHash.m $
% History:
% 001: 01-May-2011 21:52, First version.
% 007: 10-Jun-2011 10:38, [Opt.Input], binary data, complex values considered.
% 011: 26-May-2012 15:57, Fixed: Failed for binary input and empty data.
% 014: 04-Nov-2012 11:37, Consider Mex-, MDL- and P-files also.
%      Thanks to David (author 243360), who found this bug.
%      Jan Achterhold (author 267816) suggested to consider Java objects.
% 016: 01-Feb-2015 20:53, Java heap space exhausted for large files.
%      Now files are process in chunks to save memory.
% 017: 15-Feb-2015 19:40, Collsions: Same hash for different data.
%      Examples: zeros(1,1) and zeros(1,1,0)
%                complex(0) and zeros(1,1,0,0)
%      Now the number of dimensions is included, to avoid this.
% 022: 30-Mar-2015 00:04, Bugfix: Failed for strings and [] without TYPECASTX.
%      Ross found these 2 bugs, which occur when TYPECASTX is not installed.
%      If you need the base64 format padded with '=' characters, adjust
%      fBase64_enc as you like.
% 026: 29-Jun-2015 00:13, Changed hash for STRUCTs.
%      Struct arrays are analysed field by field now, which is much faster.
% 027: 13-Sep-2015 19:03, 'ascii' input as abbrev. for Input='bin' and UINT8().
% 028: 15-Oct-2015 23:11, Example values in help section updated to v022.
% 029: 16-Oct-2015 22:32, Use default options for empty input.
% 031: 28-Feb-2016 15:10, New hash value to get same reply as GetMD5.
%      New Matlab version (at least 2015b) use a fast method for TYPECAST, such
%      that calling James Tursa's TYPECASTX is not needed anymore.
%      Matlab 6.5 not supported anymore: MException for CATCH.
% 033: 18-Jun-2016 14:28, BUGFIX: Failed on empty files.
%      Thanks to Christian (AuthorID 2918599).

% OPEN BUGS:
% Nath wrote:
% function handle refering to struct containing the function will create
% infinite loop. Is there any workaround ?
% Example:
%   d= dynamicprops();
%   addprop(d,'f');
%   d.f= @(varargin) struct2cell(d);
%   DataHash(d.f) % infinite loop
% This is caught with an error message concerning the recursion limit now.

% Main function: ===============================================================
% Default options: -------------------------------------------------------------
Method    = 'MD5';
OutFormat = 'hex';
isFile    = false;
isBin     = false;

% Check number and type of inputs: ---------------------------------------------
nArg = nargin;
if nArg == 2
   if isa(Opt, 'struct') == 0   % Bad type of 2nd input:
      Error_L('BadInput2', '2nd input [Opt] must be a struct.');
   end
   
   % Specify hash algorithm:
   if isfield(Opt, 'Method')  && ~isempty(Opt.Method)   % Short-circuiting
      Method = upper(Opt.Method);
   end
   
   % Specify output format:
   if isfield(Opt, 'Format') && ~isempty(Opt.Format)    % Short-circuiting
      OutFormat = Opt.Format;
   end
   
   % Check if the Input type is specified - default: 'array':
   if isfield(Opt, 'Input') && ~isempty(Opt.Input)      % Short-circuiting
      if strcmpi(Opt.Input, 'File')
         if ischar(Data) == 0
            Error_L('CannotOpen', '1st input FileName must be a string');
         end
         isFile = true;
         
      elseif strncmpi(Opt.Input, 'bin', 3)  % Accept 'binary' also
         if (isnumeric(Data) || ischar(Data) || islogical(Data)) == 0 || ...
               issparse(Data)
            Error_L('BadDataType', ...
               '1st input must be numeric, CHAR or LOGICAL for binary input.');
         end
         isBin = true;
         
      elseif strncmpi(Opt.Input, 'asc', 3)  % 8-bit ASCII characters
         if ~ischar(Data)
            Error_L('BadDataType', ...
               '1st input must be a CHAR for the input type ASCII.');
         end
         isBin = true;
         Data  = uint8(Data);
      end
   end
   
elseif nArg == 0  % Reply version of this function:
   R = Version_L;
   
   if nargout == 0
      disp(R);
   else
      Hash = R;
   end
   
   return;
   
elseif nArg ~= 1  % Bad number of arguments:
   Error_L('BadNInput', '1 or 2 inputs required.');
end

% Create the engine: -----------------------------------------------------------
try
   Engine = java.security.MessageDigest.getInstance(Method);
catch
   Error_L('BadInput2', 'Invalid algorithm: [%s].', Method);
end

% Create the hash value: -------------------------------------------------------
if isFile
   % Open the file:
   FID = fopen(Data, 'r');
   if FID < 0
      % Check existence of file:
      Found = FileExist_L(Data);
      if Found
         Error_L('CantOpenFile', 'Cannot open file: %s.', Data);
      else
         Error_L('FileNotFound', 'File not found: %s.', Data);
      end
   end
   
   % Read file in chunks to save memory and Java heap space:
   Chunk = 1e6;      % Fastest for 1e6 on Win7/64, HDD
   Count = Chunk;    % Dummy value to satisfy WHILE condition
   while Count == Chunk
      [Data, Count] = fread(FID, Chunk, '*uint8');
      if Count ~= 0  % Avoid error for empty file
         Engine.update(Data);
      end
   end
   fclose(FID);
   
   % Calculate the hash:
   Hash = typecast(Engine.digest, 'uint8');
   
elseif isBin             % Contents of an elementary array, type tested already:
   if isempty(Data)      % Nothing to do, Engine.update fails for empty input!
      Hash = typecast(Engine.digest, 'uint8');
   else                  % Matlab's TYPECAST is less elegant:
      if isnumeric(Data)
         if isreal(Data)
            Engine.update(typecast(Data(:), 'uint8'));
         else
            Engine.update(typecast(real(Data(:)), 'uint8'));
            Engine.update(typecast(imag(Data(:)), 'uint8'));
         end
      elseif islogical(Data)               % TYPECAST cannot handle LOGICAL
         Engine.update(typecast(uint8(Data(:)), 'uint8'));
      elseif ischar(Data)                  % TYPECAST cannot handle CHAR
         Engine.update(typecast(uint16(Data(:)), 'uint8'));
         % Bugfix: Line removed
      end
      Hash = typecast(Engine.digest, 'uint8');
   end
else                 % Array with type:
   Engine = CoreHash(Data, Engine);
   Hash   = typecast(Engine.digest, 'uint8');
end

% Convert hash specific output format: -----------------------------------------
switch OutFormat
   case 'hex'
      Hash = sprintf('%.2x', double(Hash));
   case 'HEX'
      Hash = sprintf('%.2X', double(Hash));
   case 'double'
      Hash = double(reshape(Hash, 1, []));
   case 'uint8'
      Hash = reshape(Hash, 1, []);
   case 'base64'
      Hash = fBase64_enc(double(Hash));
   otherwise
      Error_L('BadOutFormat', ...
         '[Opt.Format] must be: HEX, hex, uint8, double, base64.');
end
end
% return;

% ******************************************************************************
function Engine = CoreHash(Data, Engine)
% This methods uses the slower TYPECAST of Matlab

% Consider the type and dimensions of the array to distinguish arrays with the
% same data, but different shape: [0 x 0] and [0 x 1], [1,2] and [1;2],
% DOUBLE(0) and SINGLE([0,0]):
% <  v016: [class, size, data]. BUG! 0 and zeros(1,1,0) had the same hash!
% >= v016: [class, ndims, size, data]
Engine.update([uint8(class(Data)), ...
              typecast(uint64([ndims(Data), size(Data)]), 'uint8')]);
           
if issparse(Data)                    % Sparse arrays to struct:
   [S.Index1, S.Index2, S.Value] = find(Data);
   Engine                        = CoreHash(S, Engine);
elseif isstruct(Data)                % Hash for all array elements and fields:
   F = sort(fieldnames(Data));       % Ignore order of fields
   for iField = 1:length(F)          % Loop over fields
      aField = F{iField};
      Engine.update(uint8(aField));
      for iS = 1:numel(Data)         % Loop over elements of struct array
         Engine = CoreHash(Data(iS).(aField), Engine);
      end
   end
elseif iscell(Data)                  % Get hash for all cell elements:
   for iS = 1:numel(Data)
      Engine = CoreHash(Data{iS}, Engine);
   end
elseif isempty(Data)                 % Nothing to do
elseif isnumeric(Data)
   if isreal(Data)
      Engine.update(typecast(Data(:), 'uint8'));
   else
      Engine.update(typecast(real(Data(:)), 'uint8'));
      Engine.update(typecast(imag(Data(:)), 'uint8'));
   end
elseif islogical(Data)               % TYPECAST cannot handle LOGICAL
   Engine.update(typecast(uint8(Data(:)), 'uint8'));
elseif ischar(Data)                  % TYPECAST cannot handle CHAR
   Engine.update(typecast(uint16(Data(:)), 'uint8'));
elseif isa(Data, 'function_handle')
   Engine = CoreHash(ConvertFuncHandle(Data), Engine);
elseif (isobject(Data) || isjava(Data)) && ismethod(Data, 'hashCode')
   Engine = CoreHash(char(Data.hashCode), Engine);
else  % Most likely a user-defined object:
   try
      BasicData = ConvertObject(Data);
   catch ME
      error(['JSimon:', mfilename, ':BadDataType'], ...
         '%s: Cannot create elementary array for type: %s\n  %s', ...
         mfilename, class(Data), ME.message);
   end
   
   try
      Engine = CoreHash(BasicData, Engine);
   catch ME
      if strcmpi(ME.identifier, 'MATLAB:recursionLimit')
         ME = MException(['JSimon:', mfilename, ':RecursiveType'], ...
            '%s: Cannot create hash for recursive data type: %s', ...
            mfilename, class(Data));
      end
      throw(ME);
   end
end
end
% return;

% ******************************************************************************
function FuncKey = ConvertFuncHandle(FuncH)
%   The subfunction ConvertFuncHandle converts function_handles to a struct
%   using the Matlab function FUNCTIONS. The output of this function changes
%   with the Matlab version, such that DataHash(@sin) replies different hashes
%   under Matlab 6.5 and 2009a.
%   An alternative is using the function name and name of the file for
%   function_handles, but this is not unique for nested or anonymous functions.
%   If the MATLABROOT is removed from the file's path, at least the hash of
%   Matlab's toolbox functions is (usually!) not influenced by the version.
%   Finally I'm in doubt if there is a unique method to hash function handles.
%   Please adjust the subfunction ConvertFuncHandles to your needs.

% The Matlab version influences the conversion by FUNCTIONS:
% 1. The format of the struct replied FUNCTIONS is not fixed,
% 2. The full paths of toolbox function e.g. for @mean differ.
FuncKey = functions(FuncH);

% Include modification file time and file size. Suggested by Aslak Grinsted:
if ~isempty(FuncKey.file)
    d = dir(FuncKey.file);
    if ~isempty(d)
        FuncKey.filebytes = d.bytes;
        FuncKey.filedate  = d.datenum;
    end
end

% ALTERNATIVE: Use name and path. The <matlabroot> part of the toolbox functions
% is replaced such that the hash for @mean does not depend on the Matlab
% version.
% Drawbacks: Anonymous functions, nested functions...
% funcStruct = functions(FuncH);
% funcfile   = strrep(funcStruct.file, matlabroot, '<MATLAB>');
% FuncKey    = uint8([funcStruct.function, ' ', funcfile]);

% Finally I'm afraid there is no unique method to get a hash for a function
% handle. Please adjust this conversion to your needs.
end
% return;

% ******************************************************************************
function DataBin = ConvertObject(DataObj)
% Convert a user-defined object to a binary stream. There cannot be a unique
% solution, so this part is left for the user...

try    % Perhaps a direct conversion is implemented:
   DataBin = uint8(DataObj);
   
   % Matt Raum had this excellent idea - unfortunately this function is
   % undocumented and might not be supported in te future:
   % DataBin = getByteStreamFromArray(DataObj);
   
catch  % Or perhaps this is better:
   WarnS   = warning('off', 'MATLAB:structOnObject');
   DataBin = struct(DataObj);
   warning(WarnS);
end
end
% return;

% ******************************************************************************
function Out = fBase64_enc(In)
% Encode numeric vector of UINT8 values to base64 string.
% The intention of this is to create a shorter hash than the HEX format.
% Therefore a padding with '=' characters is omitted on purpose.

Pool = [65:90, 97:122, 48:57, 43, 47];  % [0:9, a:z, A:Z, +, /]
v8   = [128; 64; 32; 16; 8; 4; 2; 1];
v6   = [32, 16, 8, 4, 2, 1];

In  = reshape(In, 1, []);
X   = rem(floor(In(ones(8, 1), :) ./ v8(:, ones(length(In), 1))), 2);
Y   = reshape([X(:); zeros(6 - rem(numel(X), 6), 1)], 6, []);
Out = char(Pool(1 + v6 * Y));
end
% return;

% ******************************************************************************
function Ex = FileExist_L(FileName)
% A more reliable version of EXIST(FileName, 'file'):
dirFile = dir(FileName);
if length(dirFile) == 1
   Ex = ~(dirFile.isdir);
else
   Ex = false;
end
end
% return;

% ******************************************************************************
function R = Version_L()
% The output differs between versions of this function. So give the user a
% chance to recognize the version:
% 1: 01-May-2011, Initial version
% 2: 15-Feb-2015, The number of dimensions is considered in addition.
%    In version 1 these variables had the same hash:
%    zeros(1,1) and zeros(1,1,0), complex(0) and zeros(1,1,0,0)
% 3: 29-Jun-2015, Struct arrays are processed field by field and not element
%    by element, because this is much faster. In consequence the hash value
%    differs, if the input contains a struct.
% 4: 28-Feb-2016 15:20, same output as GetMD5 for MD5 sums. Therefore the
%    dimensions are casted to UINT64 at first.
R.HashVersion = 4;
R.Date        = [2016, 2, 28];

R.HashMethod  = {};
try
   Provider = java.security.Security.getProviders;
   for iProvider = 1:numel(Provider)
      S     = char(Provider(iProvider).getServices);
      Index = strfind(S, 'MessageDigest.');
      for iDigest = 1:length(Index)
         Digest       = strtok(S(Index(iDigest):end));
         Digest       = strrep(Digest, 'MessageDigest.', '');
         R.HashMethod = cat(2, R.HashMethod, {Digest});
      end
   end
catch ME
   fprintf(2, '%s\n', ME.message);
   R.HashMethod = 'error';
end
end
% return;

% ******************************************************************************
function Error_L(ID, varargin)

error(['JSimon:', mfilename, ':', ID], ['*** %s: ', varargin{1}], ...
   mfilename, varargin{2:nargin - 1});
end
% return;
end

function [v, pos] = deserialize(m, pos)
% DESERIALIZE converts the output of SERIALIZE back into a matlab object.
%
% v = DESERIALIZE(m)
% [v, pos] = DESERIALIZE(m, pos)
%
% m is the series of bytes created by SERIALIZE. Integer numeric types are
% automatically converted back to doubles. The optional input/output 'pos'
% is the position to start reading from, and is returned as pointing to the
% first unused byte.
%
% If all the data is supposed to be decoded, the following should be true
% after execution.
%
% pos == numel(m)+1
%
% By Tim Hutt, 19/11/2010
%
% Updated 16/12/2010 - Fix bug with matrices.

	if ~isnumeric(m)
		error('Input must be numeric (and uint8)');
	end
	if ~strcmp(class(m), 'uint8') 
		error('Input must be uint8');
	end
	if nargin < 2
		pos = 1;
	end
	if pos > numel(m)
		error('Input too small')
	end
	
	cls = byteToClass(m(pos));
	
	switch (cls)
		case {'double', 'single', 'logical', 'char', ...
				'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'}
			[v, pos] = deserializeMatrix(m, pos);			
		case 'struct'
			[v, pos] = deserializeStruct(m, pos);			
		case 'cell'
			[v, pos] = deserializeCell(m, pos);			
		otherwise
			error('Unknown class');
    end

    
    
function [v, pos] = deserializeMatrix(m, pos)
	cls = byteToClass(m(pos));
	pos = pos + 1;
	ndms = double(m(pos));
	pos = pos + 1;
	dms = [];
	for ii = 1:ndms
		dms(ii) = double(typecast(m(pos:pos+3), 'uint32'));
		pos = pos + 4;
	end
	
	nbytes = prod(dms) * sizeof(cls);
	
	% Data.
	switch cls
		case 'char'
			v = char(m(pos:pos+nbytes-1));
		case 'logical'
			v = logical(m(pos:pos+nbytes-1));
		otherwise
			v = double(typecast(m(pos:pos+nbytes-1), cls));
	end
	
	pos = pos + nbytes;
	v = reshape(v, [dms 1 1]);
end

function sz = sizeof(cls)
	switch cls
		case {'double', 'int64', 'uint64'}
			sz = 8;
		case {'single', 'int32', 'uint32'}
			sz = 4;
		case {'int16', 'uint16'}
			sz = 2;
		case {'logical', 'char', 'int8', 'uint8'}
			sz = 1;
		otherwise
			error('Unknown class');
	end
end

function cls = byteToClass(b)
	switch b
		case 0
			cls = 'double';
		case 1
			cls = 'single';
		case 2
			cls = 'logical';
		case 3
			cls = 'char';
		case 4
			cls = 'int8';
		case 5
			cls = 'uint8';
		case 6
			cls = 'int16';
		case 7
			cls = 'uint16';
		case 8
			cls = 'int32';
		case 9
			cls = 'uint32';
		case 10
			cls = 'int64';
		case 11
			cls = 'uint64';
		case 254
			cls = 'cell';
		case 255
			cls = 'struct';
		otherwise
			error('Unknown class');
	end
end

function [v, pos] = deserializeCell(m, pos)
	pos = pos + 1; % We know it is a cell.
	
	
	ndms = double(m(pos));
	pos = pos + 1;
	dms = [];
	for ii = 1:ndms
		dms(ii) = double(typecast(m(pos:pos+3), 'uint32'));
		pos = pos + 4;
	end
	
	nels = prod(dms);
	
	v = {};
	for ii = 1:nels
		[v{ii}, pos] = deserialize(m, pos);
	end
	v = reshape(v, [dms 1 1]);
end

function [v, pos] = deserializeStruct(m, pos)
	pos = pos + 1; % We know it is a struct.
	
	% Number of field names.
	nfields = double(m(pos));
	pos = pos + 1;
	% Field names.
	for ii = 1:nfields
		fnlen = double(typecast(m(pos:pos+3), 'uint32'));
		pos = pos + 4;
		fieldNames{ii} = char(m(pos:pos+fnlen-1)).';
		pos = pos + fnlen;
	end

	% Dimensions
	ndms = double(m(pos));
	pos = pos + 1;
	dms = [];
	for ii = 1:ndms
		dms(ii) = double(typecast(m(pos:pos+3), 'uint32'));
		pos = pos + 4;
	end
	
	nels = prod(dms);
	
	v = [];
	for ii = 1:nels
		for ff = 1:nfields
			[v(ii).(fieldNames{ff}), pos] = deserialize(m, pos);
		end
	end
	v = reshape(v, [dms 1 1]);

end

    
    
end