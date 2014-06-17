% vim:tabstop=4 expandtab syntax=tex:
%
% nih.cls - class for preparing NIH (US National Institutes of Health)
% grant applications.
% 
% Copyright (C) 2003, 2004, 2005, 2006 Bruce Donald.
% 
% This work may be distributed and/or modified under the
% conditions of the LaTeX Project Public License, either version 1.3c
% of this license or (at your option) any later version.
% The latest version of this license is in
%   http://www.latex-project.org/lppl.txt
% and version 1.3 or later is part of all distributions of LaTeX
% version 2003/12/01 or later.
%
% This work has the LPPL maintenance status "author-maintained".
% 
% The Current Maintainer of this work is Bruce Donald
% (http://www.cs.duke.edu/brd/).
%
% This work consists of only this one file, nih.cls.
% Accompanying example files are not under the LPPL.
%
% last changed by Bruce Donald in 2005
% to mimic the two PHS forms.
%
% Previous modified by Alan Shields, May 2003
% email: Alan-Shields2omrf.ouhsc.edu
% Heavily based upon (read: copied almost exactly from)
%    nih.sty by Bruce R. Donald and Matthew L. Fidler
%
\NeedsTeXFormat{LaTeX2e}
%\ProvidesClass{nih}[2005 NIH PHS 398 as of June 2005]
\ProvidesClass{nih_mod}[2005/06/01 NIH PHS 398 as of June 2005]

% includes needed for preliminary code
\RequirePackage{ifthen}

% Preliminary code
\newcommand{\nih@kk}{~}

\newboolean{nih@blank}
\setboolean{nih@blank}{false}% not strictly necessary, but let's be clear about it

% Options
%%%% Changed by M A Lewis, Ph.D. (mal11 at alumni.cwru.edu)
%%%% NIH suggests Helvetica or Arial 12pt. Use Helvetica Postscript fonts for 12pt also.
%\DeclareOption{12pt}{\AtEndOfClass{\RequirePackage{helvet}\renewcommand{\familydefault}{\sfdefault}}}
%\DeclareOption{12pt}{\AtEndOfClass{\RequirePackage{times}}}
%\DeclareOption{12pt}{\AtEndOfClass{\RequirePackage{arial}}}
\DeclareOption{11pt}{\AtEndOfClass{\RequirePackage{helvet}\renewcommand{\familydefault}{\sfdefault}}}
\DeclareOption{nihblank}{\setboolean{nih@blank}{true}}
\DeclareOption{kk}{\renewcommand{\nih@kk}{\small KK}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax

% The class we're overwriting
\LoadClassWithOptions{article}


% Begin Requirements
%%%% Changed by M A Lewis, Ph.D. (mal11 at alumni.cwru.edu)
%%%% PHS 398 allows 0.5in magins everywhere. Simplify page layout using geometry package.
%\RequirePackage{fullpage}
%\RequirePackage[letterpaper,margin=0.5in,includeheadfoot]{geometry}
\RequirePackage[letterpaper,margin=0.5in]{geometry}
\RequirePackage{fancyhdr}
\RequirePackage{amsmath}



%%%% More code
% preamble stuff
\newcommand{\nih@PIname}{Vardavas, Raffaele}
\newcommand{\piname}[1]{\renewcommand{\nih@PIname}{#1}}

% constants
\newcommand{\nih@sillysize}{\scriptsize}

% set the page style
\pagestyle{empty}
%\pagestyle{fancy}
%\chead{\fancyplain{}{%
%        \ifthenelse{\boolean{nih@blank}}%
%            {}% fi
%            {%
%                {{\sf\scriptsize Principal Investigator/Program Director%
%                    (Last, First, Middle):}}%
%            }% esle
%    }
%}
%\rhead{\fancyplain{}{%
%        \ifthenelse{\boolean{nih@blank}}%
%            {\ \ \ {\nih@PIname}\ \ \ \ \ }% fi
%            {%
%                \underline{\ \ \ {\nih@PIname}\ \ \ \ \ }%
%            }% esle
%    }
%}
%\lhead{\fancyplain{}{\nih@kk}}
%\lfoot{\fancyplain{}{%
%        \ifthenelse{\boolean{nih@blank}}%
%            {}% fi
%            {\sf\footnotesize PHS 398/2590 (Rev.~09/04)\\}% esle
%    }
%}
%\cfoot{\ifthenelse{\boolean{nih@blank}}%
%        {\thepage}% fi
%        {\sf\footnotesize{Page}}% esle
%    \ \underline{\ \ \ {\sf\thepage}\ \ \ }%
%}
%\rfoot{\fancyplain{}%
%    {\sf\footnotesize{\textbf{Continuation Format Page}}}%
%}
%
%\renewcommand{\headrulewidth}{0pt}
%\renewcommand{\footrulewidth}{0.75pt}

%%%% Changed by M A Lewis, Ph.D. (mal11 at alumni.cwru.edu)
%%%% Simplify page layout by using geometry package above.
% page dimensions (note: assumes 8.5x11)
%\setlength{\oddsidemargin}{31pt}
%\setlength{\textwidth}{8.5in}
%\addtolength{\textheight}{0in} % 0.5in margins on both sides
%\addtolength{\textwidth}{-1.1in}
%\setlength{\oddsidemargin}{-1in}
%\addtolength{\oddsidemargin}{0.505in}
%\setlength{\evensidemargin}{\oddsidemargin}
%\addtolength{\textheight}{0.2in}
%\addtolength{\topmargin}{-0.2in}
%\renewcommand{\baselinestretch}{.9}
%\headwidth=\textwidth

\addtolength{\headheight}{2.5pt}

% rename the bibliography section
\AtBeginDocument{\renewcommand{\refname}{Literature~Cited}}
%FIXME: something is going on with the bibliography style. Dunno what.


\def\Paper{grant application}
\def\paper{application}
\def\refappendix{Sec.}

\def\poster{(Poster)}

%Note from brd
\long\def\todo#1{{\bf{To do:}} #1}
%\long\def\todo#1{}
\def\ICRA{IEEE International Conference on Robotics and Automation (ICRA)}

\long\def\squeezable#1{#1}

%\def\a5{$\alpha_{_5}$

\def\a5{5}

%\def\mycaptionsize{\normalsize}
%\def\mycaptionsize{\small}
%\def\mycaptionsize{\small}
\def\mycaptionsize{\footnotesize}
\def\mycodesize{\footnotesize}
\def\myeqnsize{\small}

\def\sheading#1{{\bf #1:}\ }
\def\sheading#1{\subsubsection{#1}}
%\def\sheading#1{\bigskip {\bf #1.}}

\def\ssheading#1{\noindent {\bf #1.}\ }

\newtheorem{hypothesis}{Hypothesis}
\long\def\hyp#1{\begin{hypothesis} #1 \end{hypothesis}}

\newcommand{\rv}[1]{\textbf{\color{red}{ } } {\color{red} RV comment:} \color{blue}{ #1}}
\newcommand{\sn}[1]{\textbf{\color{green}{ } } {\color{green} SN comment:} \color{blue}{ #1}}
\newcommand{\field}[1]{\mathbb{#1}}


\def\cbk#1{[{\em #1}]}

\def\R{\mathbb{R}}
\def\midv{\mathop{\,|\,}}
\def\Fscr{\mathcal{F}}
\def\Gscr{\mathcal{G}}
\def\Sscr{\mathcal{S}}
\def\set#1{{\{#1\}}}
\def\edge{\!\rightarrow\!}
\def\dedge{\!\leftrightarrow\!}
\newcommand{\EOP}{\nolinebreak[1]~~~\hspace*{\fill} $\Box$\vspace*{\parskip}\vspace*{1ex}}
%my way of doing starred references
\newcommand{\mybibitem}[1]{\bibitem{#1}
\label{mybiblabel:#1}}
\newcommand{\BC}{[}
\newcommand{\EC}{]}
\newcommand{\mycite}[1]{\ref{mybiblabel:#1}\nocite{#1}}
\newcommand{\starcite}[1]{\ref{mybiblabel:#1}\citeusmark\nocite{#1}}


\def\degree{$^\circ$}
\def\R{\mathbb{R}}
\def\Fscr{\mathcal{F}}
\def\set#1{{\{#1\}}}
\def\edge{\!\rightarrow\!}
\def\dedge{\!\leftrightarrow\!}

\long\def\gobble#1{}
\def\Jigsaw{{\sc Jigsaw}}
\def\ahelix{\ensuremath{\alpha}-helix}
\def\ahelices{\ensuremath{\alpha}-helices}
\def\ahelical{$\alpha$-helical}
\def\bstrand{\ensuremath{\beta}-strand}
\def\bstrands{\ensuremath{\beta}-strands}
\def\bsheet{\ensuremath{\beta}-sheet}
\def\bsheets{\ensuremath{\beta}-sheets}
\def\hone{{\ensuremath{^1}\rm{H}}}
\def\htwo{{$^{2}$H}}
\def\cthir{{\ensuremath{^{13}}\rm{C}}}
\def\nfif{{\ensuremath{^{15}}\rm{N}}}
\def\hn{{\rm{H}\ensuremath{^\mathrm{N}}}}
\def\hnone{{\textup{H}\ensuremath{^1_\mathrm{N}}}}
\def\ca{{\rm{C}\ensuremath{^\alpha}}}
\def\catwel{{\ensuremath{^{12}}\rm{C}\ensuremath{^\alpha}}}
\def\ha{{\rm{H}\ensuremath{^\alpha}}}
\def\cb{{\rm{C}\ensuremath{^\beta}}}
\def\hb{{\rm{H}\ensuremath{^\beta}}}
\def\hg{{\rm{H}\ensuremath{^\gamma}}}
\def\dnn{{\ensuremath{d_{\mathrm{NN}}}}}
\def\dan{{\ensuremath{d_{\alpha \mathrm{N}}}}}
\def\jconst{{\ensuremath{^{3}\mathrm{J}_{\mathrm{H}^{\mathrm{N}}\mathrm{H}^{\alpha}}}} }
\def\cbfb{{CBF-$\beta$}}

\newtheorem{defn}{Definition}
\newtheorem{claim}{Claim}

    \gobble{
    \psfrag{CO}[][]{\colorbox{white}{C}}
    \psfrag{OO}[][]{\colorbox{white}{O}}
    \psfrag{CA}[][]{\colorbox{white}{\ca}}
    \psfrag{HA}[][]{\colorbox{white}{\ha}}
    \psfrag{CB}[][]{\colorbox{white}{\cb}}
    \psfrag{HB}[][]{\colorbox{white}{\hb}}
    \psfrag{HN}[][]{\colorbox{white}{\hn}}
    \psfrag{N15}[][]{\colorbox{white}{\nfif}}
    \psfrag{dnn}[][]{\dnn}
    \psfrag{dan}[][]{\dan}
    \psfrag{phi}[][]{$\phi$}
    }

\newenvironment{closeenumerate}{\begin{list}{\arabic{enumi}.}{\topsep=0in\itemsep=0in\parsep=0in\usecounter{enumi}}}{\end{list}}
\def\CR{\hspace{0pt}}           % ``invisible'' space for line break



\newif\ifdbspacing
%\dbspacingtrue  % For double spacing
\dbspacingfalse  % For normal spacing

\ifdbspacing
 \doublespacing
 \newcommand{\capspacing}{\doublespace\mycaptionsize}
\else
 \newcommand{\capspacing}{\mycaptionsize}
\fi

\def\rulefigure{\smallskip\hrule}

% \def\codesize{\normalsize}
\def\codesize{\small}

% Can use macros \be, \ee, \en as shortcuts
%  for \begin{enumerate}, \end{enumerate}, \item
%  respectively.

\def\be{\begin{enumerate}}   % Begin Enumerate
\def\ee{\end{enumerate}}     % End Enumerate
\def\en{\item}               % ENtry (item)
\def\bi{\begin{itemize}}     % Begin Itemize
\def\ei{\end{itemize}}       % End Itemize
\def\bv{\begin{verbatim}}    % Begin Verbatim
\def\ev{\end{verbatim}}      % End Verbatim

\def\matlab{{\sc matlab} }
\def\amber{{\sc amber} }
\def\KS{{$K^*$}}
\def\KSM{{K^*}} % K-Star Math
\def\KSTM{{\tilde{K}^*}}  % K-Star Tilde Math (appx K*)
\def\KOP{{$K^{\dagger}_{o}$}}  % K-Star Optimal partial
\def\KOPM{{K^{\dagger}_{o}}}  % K-Star Optimal partial Math
\def\KP{{$K^{\dagger}$}}  % K-Star partial
\def\KPM{{K^{\dagger}}}  % K-Star partial Math
\def\KTPM{{\tilde{K}^{\dagger}}}  % K-Star Tilde partial Math
\def\KD{{$K_{_D}$}}
\def\KA{{$K_{_A}$}}
\def\qpM{{q_{_P}}}
\def\qlM{{q_{_L}}}
\def\qplM{{q_{_{PL}}}}
\def\qSplM{{q^*_{_{PL}}}}
\def\KSO{{$K^*_{o}$}} % K-Star Optimal
\def\KSOM{{K^*_{o}}}  % K-Star Optimal Math
\def\CBFB{{CBF-$\beta$}}   % Core binding factor beta
\def\argmin{\mathop{\mathrm{argmin}}}
\def\rhl#1{{\em \underline{RYAN}: *\{{#1}\}*}}
\def\set#1{{\left\{ #1 \right\}}}
\def\Escr{{\mathcal{E}}}
\def\Jscr{{\mathcal{J}}}
\def\Kscr{{\mathcal{K}}}
\def\th{{$^{{\mathrm{th}}}$}}

\newtheorem{proposition}{Proposition}
\newtheorem{lemma}{Lemma}

