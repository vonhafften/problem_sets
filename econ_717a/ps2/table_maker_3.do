
global raw_data_stddiff = round($raw_data_stddiff, 0.01)
global matched_stddiff = round($matched_stddiff, 0.01)

global table_name = "table_$table_number.tex"

local row1 = "Difference & $raw_data_diff & $matched_diff \\" 
local row2 = "SE & $att_low_se & $att_med_se \\"

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
