function [h]=boxplot(x,y,varargin)
  ## Usage: boxplot(x,y,KEY,VALUE)
  ##
  ## KEYS (long and short)
  ## ----
  ## probability
  ##          p: a vector of cumulative
  ##             probabiity values for
  ##             quantile calculation
  ##             default is
  ##             p = [0,0.25,0.5,0.75,1]
  ## quantile
  ##          q: pre-calculated
  ##             quantiles (p not needed)
  ##
  ## box width, width
  ##          w: a factor for the width
  ##             of the boxes (relative)
  ##             defaults to 1.0
  ##
  ##
  ## face color,fill,patch
  ##         fc: box fill color
  ##             defaults to white
  ##
  ## edge color, box color
  ##         bc: box edge color
  ##             defaults to blue
  ##
  ## median color, median
  ##         mc: color of the median line
  ##             defaults to red
  ##
  ## KEYS can have spaces and capital letters in them
  n=nargin-2; # (key,value) pairs
  p=[0,0.25,0.5,0.75,1]; #default
  bc="blue";
  mc="red";
  fc="white";
  BoxWidth=1;
  assert(mod(n,2)==0);
  Q=[];
   for i=1:2:n
    switch(strrep(lower(varargin{i})," ",""));
      case {"p","probability","cumulativeprobability"}
	p=varargin{i+1};
	assert(length(p)==5);
      case {"q","quantile"}
	Q=varargin{i+1}';
      case {"w","width","boxwidth"}
	BoxWidth = varargin{i+1};
	assert(isscalar(BoxWidth));
      case {"facecolor","fill","patch","fc"}
	fc=varargin{i+1};
	assert(isvector(fc) && length(fc)>=3);
      case {"edgecolor","boxcolor","bc","ec"}
	bc=varargin{i+1};
	assert(isvector(bc) && length(bc)>=3);
      case {"median","mediancolor","mc"}
	mc=varargin{i+1};
	assert(isvector(mc) && length(mc)>=3);
      otherwise
	error("unknown option: %s",varargin{i});
    endswitch
  endfor

  if isempty(Q)
    assert(isvector(x) && ismatrix(y) && length(x)==size(y,2));
    Q=quantile(y,p)';
  else
    assert(size(Q,2)==5 && size(Q,1)==length(x));
  endif
  h=struct();
  D=diff(x)*BoxWidth;
  D=cat(2,D,D(end));
  #clf;
  hold on;
  for j=1:size(Q,1)
    d=D(j)*0.4;
    h(j).whisker=plot(x(j)*ones(1,2),Q(j,[1,5]),"color",bc);
    h(j).box=patch(x(j)+[-d,+d,+d,-d],[Q(j,[2,2]),Q(j,[4,4])],"edgecolor",bc,"facecolor",fc);
    h(j).median=plot(x(j)+[-d,+d],Q(j,[3,3]),"color",mc);
  endfor
  hold off;
endfunction
