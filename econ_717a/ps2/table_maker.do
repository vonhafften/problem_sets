global unm_coef = round($unm_coef, 0.01)
global unm_se = round($unm_se, 0.01)

global att_coarse = round($att_coarse, 0.01)
global att_fine = round($att_fine, 0.01)

global att_coarse_se = round($att_coarse_se, 0.01)
global att_fine_se = round($att_fine_se, 0.01)

global table_name = "table_$table_number.tex"

local row1 = "Difference & $unm_coef & $att_coarse & $att_fine \\" 
local row2 = "SE & $unm_se & $att_coarse_se & $att_fine_se \\"

texdoc init $table_name, replace

/*tex
\begin{table}[h!]
\begin{center}
\begin{tabular}{lrrr}
\toprule
& Unmatched & ATT for \texttt{pscorea} & ATT for \texttt{pscoreb}  \\
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
