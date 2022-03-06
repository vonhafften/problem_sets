
global att_low = round($att_low, 0.01)
global att_med = round($att_med, 0.01)
global att_hi = round($att_hi, 0.01)

global att_low_se = round($att_low_se, 0.01)
global att_med_se = round($att_med_se, 0.01)
global att_hi_se = round($att_hi_se, 0.01)

global table_name = "table_$table_number.tex"

local row1 = "Difference & $att_low & $att_med & $att_hi \\" 
local row2 = "SE & $att_low_se & $att_med_se & $att_hi_se \\"

texdoc init $table_name, replace

/*tex
\begin{table}[h!]
\begin{center}
\begin{tabular}{lrrr}
\toprule
 & Bandwidth = 0.02 & Bandwidth = 0.2 & Bandwidth = 2.0  \\
\hline
tex*/ 

texdoc write `row1'
texdoc write `row2'

/*tex
\bottomrule
\end{tabular}
\end{center}
\end{table}
tex*/


texdoc close 
