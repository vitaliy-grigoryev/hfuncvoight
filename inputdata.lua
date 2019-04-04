all_a   = {0.005, 0.01, 0.05, 0.1} -- a
all_b   = {0, 1e-4, 1e-3, 1e-2} -- beta
all_l   = {1e-1, 1e-2, 1e-3, 0} -- 1-lambda
all_x0  = {0.0, 3.0, 7.0, 10.0} -- x0
all_eta = {0,    30,  45,   60} -- eta2
all_eta1= {80} -- eta1

ak   = 2
bk   = 1
lk   = 4
x0k  = 3
etak = 1
eta1 = 1

params = {"b"}
-- params= {"a", "b", "l", "x0", "eta"}
types = {"oc"}	--, "ml", "ll", "mc", "cc", "la", "oc"}
-- types = {"bl", "bc"}

angles = {} 
angles["ml"] = "\\theta_0"
angles["la"] = "\\theta_0"
angles["ll"] = "\\Delta \\theta"
angles["mc"] = "\\theta_0"
angles["oc"] = "" -- no angles
angles["cc"] = "\\Delta \\theta"

xx = {}
xx["ml"] = "x_0"
xx["la"] = "X"
xx["ll"] = "X"
xx["mc"] = "X"
xx["oc"] = "x_0"
xx["cc"] = "X"

part1 = {}
part1["ml"] = [[
\begin{tikzpicture}
\def\multiplier{1}
\begin{axis}[
xlabel={$x$},
%ylabel={$I_\l(x, \eta, x_0, \eta_0) \times	\pgfmathprintnumber[sci  precision=0]{\multiplier}$},
%ylabel={$I_\l(x, \eta, x_0, \eta_0)$},
ylabel={$r_\l(x, \eta, x_0, \eta_0)$},
grid=major,
xmin=0.0,
ymin=0,
xmax=10,
%ymax=0.5,
width=8cm,
height=8cm,
restrict y to domain= 0:500,
restrict x to domain= 0:150,
]]

part1["la"] = [[
\begin{tikzpicture}
\def\multiplier{1}
\begin{axis}[
xlabel={$x$},
%ylabel={$I_\l(x, \eta, [0;X], \eta_0) \times	\pgfmathprintnumber[sci  precision=0]{\multiplier}$},
%ylabel={$\ovl {I}_\l(x, \eta)$},
ylabel={$\ovl{r}_\l(x, \eta, X, \eta_0)$},
grid=major,
xmin=0.0,
ymin=0,
xmax=10,
%ymax=0.5,
width=8cm,
height=8cm,
restrict y to domain= 0:250,
restrict x to domain= 0:150,
]]

part1["ll"] = [[
\begin{tikzpicture}
\def\multiplier{1}
\begin{axis}[
xlabel={$x$},
%ylabel={$I_\l(x, \eta, x_0, \eta_0) \times	\pgfmathprintnumber[sci  precision=0]{\multiplier}$},
%ylabel={$\ovl {I}_\l(x, \eta)$},
ylabel={$\ovl{r}_\l(x, \eta)$},
grid=major,
xmin=0.0,
ymin=0,
xmax=10,
%ymax=0.5,
width=8cm,
height=8cm,
restrict y to domain= 0:250,
restrict x to domain= 0:150,
]]

part1["cc"] = [[
\begin{tikzpicture}
\def\multiplier{1}
\begin{axis}[
xlabel={$x$},
%ylabel={$I_\cc(x, \eta, x_0, \eta_0) \times	\pgfmathprintnumber[sci  precision=0]{\multiplier}$},
%ylabel={$\ovl {I}_\cc(x, \eta)$},
ylabel={$\ovl{r}_\cc(x, \eta)$},
grid=major,
xmin=0.0,
ymin=0,
xmax=10,
%ymax=0.5,
width=8cm,
height=8cm,
restrict y to domain= 0:250,
restrict x to domain= 0:150,
]]

part1["oc"] = [[
\begin{tikzpicture}
\def\multiplier{1}
\begin{axis}[
xlabel={$x$},
%ylabel={$I_\cc(x, \eta, x_0, \eta_0) \times	\pgfmathprintnumber[sci  precision=0]{\multiplier}$},
%ylabel={$\ovl {I}_\cc(x, \eta)$},
ylabel={$r_\cc(x, \eta)$},
grid=major,
xmin=0.0,
ymin=0,
xmax=10,
%ymax=0.5,
width=8cm,
height=8cm,
restrict y to domain= 0:250,
restrict x to domain= 0:150,
]]

part1["mc"] = [[
\begin{tikzpicture}
\def\multiplier{1}
\begin{axis}[
xlabel={$x$},
%ylabel={$I_\cc(x, \eta, x_0, \eta_0) \times	\pgfmathprintnumber[sci  precision=0]{\multiplier}$},
%ylabel={$I_\cc(x, \eta, [0;X], \eta_0)$},
ylabel={$r_\cc(x, \eta, [0;X], \eta_0)$},
grid=major,
xmin=0.0,
ymin=0,
xmax=10,
%ymax=0.5,
width=8cm,
height=8cm,
restrict y to domain= 0:10050,
restrict x to domain= 0:150,
]]

part2 = [[
]

%\node[anchor=west] at (3,120) {$0.001$};
%\node[anchor=west] at (5,40) {$0.01$}; %[anchor=south]
%\node at (3.5,15) {$0.1$};
%\node at (3.5,5) {$1.0$};
%\draw (12,300) -- (9,272);
]]

part3 = [[
\pgfplotstablecreatecol[create col/expr={\thisrowno{2}*\multiplier}]{new1}{\dataa}
\pgfplotstablecreatecol[create col/expr={\thisrowno{2}*\multiplier}]{new1}{\datab}
\pgfplotstablecreatecol[create col/expr={\thisrowno{2}*\multiplier}]{new1}{\datac}
\pgfplotstablecreatecol[create col/expr={\thisrowno{2}*\multiplier}]{new1}{\datad}


\addplot[no markers, red, very thick] table[x index=0, y=new1]  {\dataa}; 
\addplot[no markers, blue!50!black, very thick] table[x index=0, y=new1]  {\datab}; %dashed
\addplot[no markers, green!50!black, very thick] table[x index=0, y=new1]  {\datac}; %dashdotted
\addplot[no markers, purple!50!black, very thick] table[x index=0, y=new1]  {\datad}; %dashdashdotted
\end{axis}
\end{tikzpicture}

]]
datas = {"\\dataa", "\\datab", "\\datac", "\\datad"}


standart = {
[0] = 	function(typ) 
			stra  = "10^{"..math.log(all_a[ak],10).."}"
			ai = ak
			strb  = "10^{"..math.log(all_b[bk],10).."}"
			bi = bk
			strl  = "10^{"..math.log(all_l[lk],10).."}"
			li = lk
			strx0 = " "..all_x0[x0k].." "
			x0i = x0k
			if typ == "ml" or typ == "mc" or typ == "la" then
				streta= " "..all_eta[etak].." "
			else 
				if typ == "oc" then
					streta= " "
				else
					streta= " "..all_eta1[eta1] - all_eta[etak].." "
				end -- if
			end -- if
			etai = etak
		end,
["a"] = function(i)
			standart[0](i)
			stra = "..."
			ai = i
		end,
["b"] = function(i)
			standart[0](i)
			strb = "..."
			bi = i
		end,
["l"] = function(i)
			standart[0](i)
			strl = "..."
			li = i
		end,
["x0"] = function(i)
			standart[0](i)
			strx0 = "..."
			x0i = i
		end,
["eta"] =function(i)
			standart[0](i)
			streta = "..."
			etai = i
		end
}

for ityp, typ in ipairs(types) do
	for ipar, par in ipairs(params) do
		print("-------->>> LUA: Different "..par)
		print("-------->>> LUA: Calculating type ["..typ.."]")
		t = io.open("out/diff_"..typ.."_"..par..".tex","w")
		standart[par](typ)
		
		t:write(part1[typ])
		t:write('title={$a = '..stra..'$, $\\beta='..strb..'$, $1-\\lambda='..strl..'$, $'..xx[typ]..'='..strx0..'$, $'..angles[typ].."="..streta..'\\deg$},')
		t:write(part2)
		for i = 1, 4 do 
			f = io.open("in.txt","w")
			standart[par](i)
			f:write(all_a[ai],  '\n')
			f:write(all_b[bi],  '\n')
			f:write(all_l[li],  '\n')
			f:write(all_x0[x0i],'\n')
			f:write(all_eta1[eta1],'\n')
			f:write(all_eta[etai],'\n')
			f:write(typ,'\n')
			f:flush()
			f:close()
	
			q = os.execute('./run')
	
			tt = io.open("files_created.txt", "r")
			t:write("\\pgfplotstableread{"..string.sub(tt:read("*l"),5).."}{"..datas[i].."}\n")
			tt:close()
		end
	
		t:write(part3)
		t:flush()
		t:close()
		-- print("-------->>> LUA: Type "..typ.." done!")	
		print("-------->>> LUA: Different "..par.." done!")	
	end -- typ
	
end -- par

