<TeXmacs|1.99.6>

<style|<tuple|generic|old-spacing|old-spacing|old-spacing>>

<\body>
  <doc-data|<doc-title|Hamiltonian of a 2D quantum
  well>|<doc-author|<\author-data|<author-name|Yunhai li>>
    \;
  </author-data>>>

  <section|Hamiltonian>

  <subsection|Kinetic terms>

  According to paper (PhysRevB.32.1043), the kinetic terms take the following
  form

  <\eqnarray>
    <tformat|<table|<row|<cell|T<rsub|e>>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2m<rsup|\<ast\>><rsub|e\<perp\>>><frac|\<partial\><rsup|2>|\<partial\>z<rsup|2><rsub|e>>,>>|<row|<cell|T<rsub|h>>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2m<rsup|\<ast\>><rsub|h\<perp\>>><frac|\<partial\><rsup|2>|\<partial\>z<rsup|2><rsub|h>>,>>|<row|<cell|T<rsub|r>>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2\<mu\>><around*|[|<frac|1|r><frac|\<partial\>|\<partial\>r>+<frac|\<partial\><rsup|2>|\<partial\>r<rsup|2>>+<frac|1|r<rsup|2>><frac|\<partial\><rsup|2>|\<partial\>\<theta\><rsup|2>>|]>\<nocomma\>,>>|<row|<cell|\<mu\>>|<cell|=>|<cell|<frac|m<rsup|\<ast\>><rsub|e<big|parallel>>m<rsup|\<ast\>><rsub|h<big|parallel>>|m<rsup|\<ast\>><rsub|e<big|parallel>>+m<rsup|\<ast\>><rsub|h<big|parallel>>>.<eq-number>>>>>
  </eqnarray>

  <subsection|Potential terms>

  <\eqnarray>
    <tformat|<table|<row|<cell|V<rsub|e>>|<cell|=>|<cell|<frac|e<rsup|2>|2\<epsilon\><rsub|1>><big|sum><rsub|n\<neq\>0><frac|q<rsub|n>|<around*|\||z<rsub|e>-z<rsub|e
    n>|\|>>,>>|<row|<cell|V<rsub|h>>|<cell|=>|<cell|<frac|e<rsup|2>|2\<epsilon\><rsub|1>><big|sum><rsub|n\<neq\>0><frac|q<rsub|n>|<around*|\||z<rsub|h>-z<rsub|h
    n>|\|>>,>>|<row|<cell|V<rsub|r>>|<cell|=>|<cell|-<frac|e<rsup|2>|\<epsilon\><rsub|1>><big|sum><rsup|+\<infty\>><rsub|n=-\<infty\>><frac|q<rsub|n>|<sqrt|<around*|(|z<rsub|e>-z<rsub|n
    h>|)><rsup|2>+r<rsup|2>>>.<eq-number>>>>>
  </eqnarray>

  All the terms are taken from Nagaosa's paper.

  <section|Trial wave function>

  The trial wave function takes the form

  <\equation>
    \<Psi\><around*|(|r,z<rsub|e>,z<rsub|h>|)>=<frac|1|l
    a><sqrt|<frac|2|\<pi\>>> \<cdot\>cos<around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|e>|)>\<cdot\>cos<around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|h>|)>\<cdot\>exp<around*|(|-<frac|r|a>|)>.
  </equation>

  It can be proved that this trial wave function is normalized

  <\eqnarray>
    <tformat|<table|<row|<cell| >|<cell|>|<cell|<big|int><rsup|l><rsub|-l><big|int><rsup|l><rsub|-l><big|int><rsup|2\<pi\>><rsub|0><big|int><rsup|+\<infty\>><rsub|0>\<Psi\><rsup|\<ast\>>\<Psi\>
    r\<mathd\>r \<mathd\>\<theta\> \<mathd\>z<rsub|e>
    \<mathd\>z<rsub|h>>>|<row|<cell|>|<cell|=>|<cell|<frac|2|l<rsup|2>a<rsup|2>\<pi\>>\<cdot\><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|e>|)>
    \<mathd\>z<rsub|e>\<cdot\><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|h>|)>
    \<mathd\>z<rsub|h>\<cdot\><big|int><rsup|2\<pi\>><rsub|0>\<mathd\>\<theta\>\<cdot\><big|int><rsup|+\<infty\>><rsub|0>r
    exp<around*|(|-<frac|2r|a>|)>\<mathd\>r>>|<row|<cell|>|<cell|=>|<cell|<frac|2|l<rsup|2>a<rsup|2>\<pi\>>\<cdot\>l\<cdot\>l\<cdot\>2\<pi\>\<cdot\><frac|a<rsup|2>|4>>>|<row|<cell|>|<cell|=>|<cell|<frac|2|l<rsup|2>a<rsup|2>\<pi\>>\<cdot\><frac|l<rsup|2>a<rsup|2>\<pi\>|2>>>|<row|<cell|>|<cell|=>|<cell|1.<eq-number>>>>>
  </eqnarray>

  <section|Matrix elements>

  <subsection|Kinetic terms>

  The matrix element of <math|T<rsub|e>> is

  <\eqnarray>
    <tformat|<table|<row|<cell|\<langle\>\<Psi\><around*|\||T<rsub|e>|\|>\<Psi\>\<rangle\>>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2m<rsup|\<ast\>><rsub|e\<perp\>>>\<langle\>\<Psi\><around*|\||<frac|\<partial\><rsup|2>|\<partial\>z<rsup|2><rsub|e>>|\|>\<Psi\>\<rangle\>>>|<row|<cell|>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2m<rsup|\<ast\>><rsub|e\<perp\>>>\<cdot\>-<frac|\<pi\>|2a<rsup|2>l<rsup|4>>\<cdot\><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|e>|)>
    \<mathd\>z<rsub|e>\<cdot\><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|h>|)>
    \<mathd\>z<rsub|h>>>|<row|<cell|>|<cell|\<cdot\>>|<cell|<big|int><rsup|2\<pi\>><rsub|0>\<mathd\>\<theta\>\<cdot\><big|int><rsup|+\<infty\>><rsub|0>r
    exp<around*|(|-<frac|2r|a>|)>\<mathd\>r>>|<row|<cell|>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2m<rsup|\<ast\>><rsub|e\<perp\>>>\<cdot\>-<frac|\<pi\>|2a<rsup|2>l<rsup|4>>\<cdot\><frac|l<rsup|2>a<rsup|2>\<pi\>|2>>>|<row|<cell|>|<cell|=>|<cell|<frac|\<hbar\><rsup|2>|2m<rsup|\<ast\>><rsub|e\<perp\>>>\<cdot\><frac|\<pi\><rsup|2>|4l<rsup|2>>.<eq-number>>>>>
  </eqnarray>

  Analogously

  <\equation>
    \<langle\>\<Psi\><around*|\||T<rsub|h>|\|>\<Psi\>\<rangle\>=<frac|\<hbar\><rsup|2>|2m<rsup|\<ast\>><rsub|h\<perp\>>>\<cdot\><frac|\<pi\><rsup|2>|4l<rsup|2>>.
  </equation>

  Now we deal with <math|T<rsub|r>>. As the trial wave function is
  independent on <math|\<theta\>>, <math|T<rsub|r>> is reduced to
  <math|T<rsub|r>=-<frac|\<hbar\><rsup|2>|2\<mu\>><around*|[|<frac|1|r><frac|\<partial\>|\<partial\>r>+<frac|\<partial\><rsup|2>|\<partial\>r<rsup|2>>|]>>.

  <\eqnarray>
    <tformat|<table|<row|<cell|\<langle\>\<Psi\><around*|\||T<rsub|r>|\|>\<Psi\>\<rangle\>>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2\<mu\>>\<langle\>\<Psi\><around*|\||<frac|1|r><frac|\<partial\>|\<partial\>r>+<frac|\<partial\><rsup|2>|\<partial\>r<rsup|2>>|\|>\<Psi\>\<rangle\>>>|<row|<cell|>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2\<mu\>><around*|[|\<langle\>\<Psi\><around*|\||<frac|1|r><frac|\<partial\>|\<partial\>r>|\|>\<Psi\>\<rangle\>+\<langle\>\<Psi\><around*|\||<frac|\<partial\><rsup|2>|\<partial\>r<rsup|2>>|\|>\<Psi\>\<rangle\>|]>>>|<row|<cell|>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2\<mu\>><around*|[|\<langle\>\<Psi\><around*|\||<frac|1|r><frac|\<partial\>|\<partial\>r>|\|>\<Psi\>\<rangle\>+<frac|2|\<pi\>a<rsup|4>l<rsup|2>>\<cdot\><frac|l<rsup|2>a<rsup|2>\<pi\>|2>|]>>>|<row|<cell|>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2\<mu\>><around*|[|\<langle\>\<Psi\><around*|\||<frac|1|r><frac|\<partial\>|\<partial\>r>|\|>\<Psi\>\<rangle\>+<frac|1|a<rsup|2>>|]>.<eq-number>>>>>
  </eqnarray>

  <\eqnarray>
    <tformat|<table|<row|<cell|\<langle\>\<Psi\><around*|\||<frac|1|r><frac|\<partial\>|\<partial\>r>|\|>\<Psi\>\<rangle\>>|<cell|=>|<cell|-<frac|2|\<pi\>a<rsup|3>l<rsup|2>>*\<cdot\><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|e>|)>
    \<mathd\>z<rsub|e>\<cdot\><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|h>|)>
    \<mathd\>z<rsub|h>>>|<row|<cell|>|<cell|\<cdot\>>|<cell|<big|int><rsup|2\<pi\>><rsub|0>\<mathd\>\<theta\>\<cdot\><big|int><rsup|+\<infty\>><rsub|0>exp<around*|(|-<frac|2r|a>|)>\<mathd\>r>>|<row|<cell|>|<cell|=>|<cell|-<frac|2|\<pi\>a<rsup|3>l<rsup|2>>\<cdot\>l<rsup|2>*\<cdot\>2\<pi\>\<cdot\><frac|a|2>>>|<row|<cell|>|<cell|=>|<cell|-<frac|2|a<rsup|2>>.<eq-number>>>>>
  </eqnarray>

  <\eqnarray>
    <tformat|<table|<row|<cell|\<langle\>\<Psi\><around*|\||T<rsub|r>|\|>\<Psi\>\<rangle\>>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2\<mu\>><around*|[|\<langle\>\<Psi\><around*|\||<frac|1|r><frac|\<partial\>|\<partial\>r>|\|>\<Psi\>\<rangle\>+<frac|1|a<rsup|2>>|]>>>|<row|<cell|>|<cell|=>|<cell|-<frac|\<hbar\><rsup|2>|2\<mu\>><around*|[|-<frac|2|a<rsup|2>>+<frac|1|a<rsup|2>>|]>>>|<row|<cell|>|<cell|=>|<cell|<frac|\<hbar\><rsup|2>|2\<mu\>>\<cdot\><frac|1|a<rsup|2>>.<eq-number>>>>>
  </eqnarray>

  <subsection|Potential terms>

  The matrix element of <math|V<rsub|e>> is

  <\eqnarray>
    <tformat|<table|<row|<cell|\<langle\>\<Psi\><around*|\||V<rsub|e>|\|>\<Psi\>\<rangle\>>|<cell|=>|<cell|<frac|e<rsup|2>|2\<epsilon\><rsub|1>><big|sum><rsub|n\<neq\>0>\<langle\>\<Psi\><around*|\||<frac|q<rsub|n>|<around*|\||z<rsub|e>-z<rsub|e
    n>|\|>>|\|>\<Psi\>\<rangle\>>>|<row|<cell|>|<cell|=>|<cell|<frac|e<rsup|2>|2\<epsilon\><rsub|1>><around*|[|<big|sum><rsup|><rsub|n\<neq\>0,even>\<langle\>\<Psi\><around*|\||<frac|q<rsub|n>|<around*|\||z<rsub|e>-<around*|(|z<rsub|e>+2n
    l|)>|\|>>|\|>\<Psi\>\<rangle\>+<big|sum><rsup|><rsub|n\<neq\>0,odd>\<langle\>\<Psi\><around*|\||<frac|q<rsub|n>|<around*|\||z<rsub|e>-<around*|(|-z<rsub|e>+2n
    l|)>|\|>>|\|>\<Psi\>\<rangle\>|]>>>|<row|<cell|>|<cell|=>|<cell|<frac|e<rsup|2>|2\<epsilon\><rsub|1>><around*|[|<big|sum><rsub|n\<neq\>0,even><frac|q<rsub|n>|2<around*|\||n
    l|\|> >\<langle\>\<Psi\>\|\<Psi\>\<rangle\>+<big|sum><rsub|n\<neq\>0,odd>\<langle\>\<Psi\><around*|\||<frac|q<rsub|n>|2<around*|\||z<rsub|e>-n
    l|\|>>|\|>\<Psi\>\<rangle\>|]>.<eq-number>>>>>
  </eqnarray>

  Let's deal with the first term, which can be devided according to <math|n>

  <\eqnarray>
    <tformat|<table|<row|<cell|<big|sum><rsub|n\<neq\>0,even><frac|q<rsub|n>|2<around*|\||n
    l|\|> >\<langle\>\<Psi\>\|\<Psi\>\<rangle\>>|<cell|=>|<cell|<big|sum><rsub|n\<neq\>0,even><frac|1|2<around*|\||n
    l|\|> >\<cdot\><around*|(|<frac|\<epsilon\><rsub|1>-\<epsilon\><rsub|2>|\<epsilon\><rsub|1>+\<epsilon\><rsub|2>>|)><rsup|<around*|\||n|\|>>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|n\<less\>0,even><frac|1|-2n
    l>\<cdot\><around*|(|<frac|\<epsilon\><rsub|1>-\<epsilon\><rsub|2>|\<epsilon\><rsub|1>+\<epsilon\><rsub|2>>|)><rsup|-n>+<big|sum><rsub|n\<gtr\>0,even><frac|1|2n
    l>\<cdot\><around*|(|<frac|\<epsilon\><rsub|1>-\<epsilon\><rsub|2>|\<epsilon\><rsub|1>+\<epsilon\><rsub|2>>|)><rsup|n>>>|<row|<cell|>|<cell|=>|<cell|2<big|sum><rsub|n\<gtr\>0,even><frac|1|2n
    l>\<cdot\><around*|(|<frac|\<epsilon\><rsub|1>-\<epsilon\><rsub|2>|\<epsilon\><rsub|1>+\<epsilon\><rsub|2>>|)><rsup|n>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|n\<gtr\>0,even><frac|1|n
    l>\<cdot\><around*|(|<frac|\<epsilon\><rsub|1>-\<epsilon\><rsub|2>|\<epsilon\><rsub|1>+\<epsilon\><rsub|2>>|)><rsup|n>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|k\<gtr\>0><frac|1|2k
    l><around*|(|<frac|\<epsilon\><rsub|1>-\<epsilon\><rsub|2>|\<epsilon\><rsub|1>+\<epsilon\><rsub|2>>|)><rsup|2k>>>|<row|<cell|>|<cell|=>|<cell|-<frac|1|2l>ln<around*|(|1-q<rsup|2>|)>\<nocomma\>>>|<row|<cell|q>|<cell|=>|<cell|<frac|\<epsilon\><rsub|1>-\<epsilon\><rsub|2>|\<epsilon\><rsub|1>+\<epsilon\><rsub|2>>.<eq-number>>>>>
  </eqnarray>

  The the second term

  <\eqnarray>
    <tformat|<table|<row|<cell|<big|sum><rsub|n\<neq\>0,odd>\<langle\>\<Psi\><around*|\||<frac|q<rsub|n>|2<around*|\||z<rsub|e>-n
    l|\|>>|\|>\<Psi\>\<rangle\>>|<cell|=>|<cell|<big|sum><rsub|n\<neq\>0,odd><frac|q<rsub|n>|2>\<langle\>\<Psi\><around*|\||<frac|1|<around*|\||z<rsub|e>-n
    l|\|>>|\|>\<Psi\>\<rangle\>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|n=0,odd><frac|q<rsub|n>|2>\<cdot\><frac|2|l<rsup|2>a<rsup|2>\<pi\>>\<cdot\><frac|a<rsup|2>l\<pi\>|2><big|int><rsup|l><rsub|-l><frac|cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|e>|)>|<around*|\||z<rsub|e>-n
    l|\|>> \<mathd\>z<rsub|e>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|n\<neq\>0,odd><frac|q<rsub|n>|2l><big|int><rsup|l><rsub|-l><frac|cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|e>|)>|<around*|\||z<rsub|e>-n
    l|\|>>\<mathd\>z<rsub|e>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|n\<gtr\>0,odd><frac|q<rsub|n>|2l><big|int><rsup|l><rsub|-l><frac|cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|e>|)>|n
    l-z<rsub|e>>\<mathd\>z<rsub|e>>>|<row|<cell|>|<cell|+>|<cell|<big|sum><rsub|n\<less\>0,odd><frac|q<rsub|n>|2l><big|int><rsup|l><rsub|-l><frac|cos<rsup|2><around*|(|<frac|\<pi\>|2l>\<cdot\>z<rsub|e>|)>|z<rsub|e>-n
    l>\<mathd\>z<rsub|e>.<eq-number>>>>>
  </eqnarray>

  Unfortunately, the results of both the two terms involve transcendental
  function <em|cosine integral>. However, at least we are sure that it is not
  dependent on the parameter <math|a>. Thus
  <math|\<langle\>\<Psi\><around*|\||V<rsub|e>|\|>\<Psi\>\<rangle\>> can be
  neglected in the variation. Similar conclusions can be drawn for
  <math|\<langle\>\<Psi\><around*|\||V<rsub|h>|\|>\<Psi\>\<rangle\>>.

  The matrix element of <math|V<rsub|r>> is

  <\eqnarray>
    <tformat|<table|<row|<cell|\<langle\>\<Psi\><around*|\||V<rsub|r>|\|>\<Psi\>\<rangle\>>|<cell|=>|<cell|-<frac|e<rsup|2>|\<epsilon\><rsub|1>><big|sum><rsup|+\<infty\>><rsub|n=-\<infty\>>\<langle\>\<Psi\><around*|\||<frac|q<rsub|n>|<sqrt|<around*|(|z<rsub|e>-z<rsub|n
    h>|)><rsup|2>+r<rsup|2>>>|\|>\<Psi\>\<rangle\>>>|<row|<cell|>|<cell|=>|<cell|-<frac|e<rsup|2>|\<epsilon\><rsub|1>><big|sum><rsup|+\<infty\>><rsub|n=-\<infty\>>q<rsub|n>\<langle\>\<Psi\><around*|\||<frac|1|<sqrt|<around*|(|z<rsub|e>-z<rsub|n
    h>|)><rsup|2>+r<rsup|2>>>|\|>\<Psi\>\<rangle\>.>>>>
  </eqnarray>

  \;

  <\eqnarray>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<langle\>\<Psi\><around*|\||<frac|1|<sqrt|<around*|(|z<rsub|e>-z<rsub|n
    h>|)><rsup|2>+r<rsup|2>>>|\|>\<Psi\>\<rangle\>>>|<row|<cell|>|<cell|=>|<cell|<frac|2|l<rsup|2>a<rsup|2>\<pi\>>\<cdot\><big|int><rsup|2\<pi\>><rsub|0>\<mathd\>\<theta\>\<cdot\><big|int><rsup|l><rsub|-l><big|int><rsup|l><rsub|-l><big|int><rsup|\<infty\>><rsub|0>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|e>|2l>|)>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|h>|2l>|)><frac|r\<cdot\>exp<around*|(|-<frac|2r|a>|)>|<sqrt|<around*|(|z<rsub|e>-z<rsub|n
    h>|)><rsup|2>+r<rsup|2>>>\<mathd\>r\<mathd\>z<rsub|e>\<mathd\>z<rsub|h>>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|=>|<cell|<frac|4|l<rsup|2>a<rsup|2>><big|int><rsup|l><rsub|-l><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|e>|2l>|)>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|h>|2l>|)>\<mathd\>z<rsub|e>\<mathd\>z<rsub|h><big|int><rsup|\<infty\>><rsub|0><frac|r\<cdot\>exp<around*|(|-<frac|2r|a>|)>|<sqrt|<around*|(|z<rsub|e>-z<rsub|n
    h>|)><rsup|2>+r<rsup|2>>>\<mathd\>r>>|<row|<cell|>|<cell|=>|<cell|<frac|2|l<rsup|2>a><big|int><rsup|l><rsub|-l><big|int><rsup|l><rsub|-l><around*|[|cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|e>|2l>|)>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|h>|2l>|)>\<mathd\>z<rsub|e>\<mathd\>z<rsub|h>\<cdot\><frac|2|a><big|int><rsup|\<infty\>><rsub|0><frac|r\<cdot\>exp<around*|(|-<frac|2r|a>|)>|<sqrt|<around*|(|z<rsub|e>-z<rsub|n
    h>|)><rsup|2>+r<rsup|2>>>\<mathd\>r|]>>>|<row|<cell|>|<cell|=>|<cell|<frac|2|l<rsup|2>a><big|int><rsup|l><rsub|-l><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|e>|2l>|)>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|h>|2l>|)>G<around*|(|z<rsub|e>-z<rsub|n
    h>|)>\<mathd\>z<rsub|e>\<mathd\>z<rsub|h>,>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|G<around*|(|\<gamma\>|)>>|<cell|=>|<cell|<frac|2<around*|\||\<gamma\>|\|>|a><around*|{|<frac|\<pi\>|2><around*|[|H<rsub|1><around*|(|<frac|2<around*|\||\<gamma\>|\|>|a>|)>-N<rsub|1><around*|(|<frac|2<around*|\||\<gamma\>|\|>|a>|)>|]>-1|}>.<eq-number>>>>>
  </eqnarray>

  Finally we have

  <\eqnarray>
    <tformat|<table|<row|<cell|\<langle\>\<Psi\><around*|\||V<rsub|r>|\|>\<Psi\>\<rangle\>>|<cell|=>|<cell|-<frac|e<rsup|2>|\<epsilon\><rsub|1>><big|sum><rsup|+\<infty\>><rsub|n=-\<infty\>>q<rsub|n><around*|[|<frac|2|l<rsup|2>a><big|int><rsup|l><rsub|-l><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|e>|2l>|)>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|h>|2l>|)>G<around*|(|z<rsub|e>-z<rsub|n
    h>|)>\<mathd\>z<rsub|e>\<mathd\>z<rsub|h>|]>>>|<row|<cell|>|<cell|=>|<cell|-<frac|e<rsup|2>|\<epsilon\><rsub|1>>\<cdot\><frac|2|l<rsup|2>a><big|int><rsup|l><rsub|-l><big|int><rsup|l><rsub|-l>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|e>|2l>|)>cos<rsup|2><around*|(|<frac|\<pi\>z<rsub|h>|2l>|)><around*|[|<big|sum><rsup|+\<infty\>><rsub|n=-\<infty\>>q<rsub|n>G<around*|(|z<rsub|e>-z<rsub|n
    h>|)>|]>\<mathd\>z<rsub|e>\<mathd\>z<rsub|h>.<eq-number>>>>>
  </eqnarray>

  <subsection|One barrier case>

  If the quantum well has only one barrier layer on either side, the
  situation will be complicated. We generalize this problem in to a quantum
  well with two barriers with different dielectric constants on both sides,
  and denote the dielectric constants of the bottom barrier, quantum well and
  top barrier as <math|\<epsilon\><rsub|1>>, <math|\<epsilon\><rsub|2>> and
  <math|\<epsilon\><rsub|3>>, respectively. Also, we shift the origin of axes
  to <math|l> in order to utilize the formulae in the literature. Then the
  trial wave function becomes

  <\equation>
    \<Psi\><rprime|'><around*|(|r,z<rsub|e>,z<rsub|h>|)>=<frac|1|l
    a><sqrt|<frac|2|\<pi\>>>\<cdot\>sin<around*|(|<frac|\<pi\>z<rsub|e>|2l>|)>\<cdot\>sin<around*|(|<frac|\<pi\>z<rsub|h>|2l>|)>\<cdot\>exp<around*|(|-<frac|r|a>|)>,z<rsub|e>,z<rsub|h>\<in\><around*|[|-2l,0|]>.
  </equation>

  Then the matrix element of <math|V<rsub|r>> becomes

  <\eqnarray>
    <tformat|<table|<row|<cell|\<langle\>\<Psi\><rprime|'><around*|\||V<rsub|r>|\|>\<Psi\><rprime|'>\<rangle\>>|<cell|=>|<cell|\<langle\>\<Psi\><rprime|'><around*|\||-<frac|e<rsup|2>|\<epsilon\><rsub|2>><around*|{|\<ast\>|}>|\|>\<Psi\><rprime|'>\<rangle\>>>|<row|<cell|>|<cell|=>|<cell|-<frac|e<rsup|2>|\<epsilon\><rsub|2>>\<langle\>\<Psi\><rprime|'><around*|\||<around*|{|\<ast\>|}>|\|>\<Psi\><rprime|'>\<rangle\>,>>>>
  </eqnarray>

  <\eqnarray>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<langle\>\<Psi\><rprime|'><around*|\||<around*|{|\<ast\>|}>|\|>\<Psi\><rprime|'>\<rangle\>>>|<row|<cell|>|<cell|=>|<cell|<frac|4|l<rsup|2>a<rsup|2>><big|int><rsup|0><rsub|-2l><big|int><rsup|0><rsub|-2l>sin<rsup|2><around*|(|<frac|\<pi\>z<rsub|e>|2l>|)>sin<rsup|2><around*|(|<frac|\<pi\>z<rsub|h>|2l>|)>\<mathd\>z<rsub|e>\<mathd\>z<rsub|h><big|int><rsup|+\<infty\>><rsub|0>r
    exp<around*|(|-<frac|2r|a>|)><around*|{|\<ast\>|}>\<mathd\>r>>|<row|<cell|>|<cell|=>|<cell|<frac|2|l<rsup|2>a><big|int><rsup|0><rsub|-2l><big|int><rsup|0><rsub|-2l>sin<rsup|2><around*|(|<frac|\<pi\>z<rsub|e>|2l>|)>sin<rsup|2><around*|(|<frac|\<pi\>z<rsub|h>|2l>|)>\<mathd\>z<rsub|e>\<mathd\>z<rsub|h>\<cdot\><frac|2|a><big|int><rsup|+\<infty\>><rsub|0>r
    exp<around*|(|-<frac|2r|a>|)><around*|{|\<ast\>|}>\<mathd\>r>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|<frac|2|a><big|int><rsup|+\<infty\>><rsub|0>r
    exp<around*|(|-<frac|2r|a>|)><around*|{|\<ast\>|}>\<mathd\>r>>|<row|<cell|>|<cell|=>|<cell|<frac|2|a><big|int><rsup|\<infty\>><rsub|0>r
    exp<around*|(|-<frac|2r|a>|)><around*|{|<frac|1|R<around*|(|z<rsub|h>|)>>+<frac|k|R<around*|(|-z<rsub|h>|)>>+<big|sum><rsup|+\<infty\>><rsub|n=1>\<lambda\><rsup|n><around*|[|<frac|k|R<around*|(|z<rsub|+n>|)>>+<frac|1|R<around*|(|-z<rsub|+n>|)>>+<frac|1|R<around*|(|z<rsub|-n>|)>>+<frac|1|k
    R<around*|(|-z<rsub|-n>|)>>|]>|}>\<mathd\>r>>|<row|<cell|>|<cell|=>|<cell|G<around*|(|z<rsub|e>-z<rsub|h>|)>+k\<cdot\>G<around*|(|z<rsub|e>+z<rsub|h>|)>>>|<row|<cell|>|<cell|+>|<cell|<big|sum><rsup|+\<infty\>><rsub|n=1>\<lambda\><rsup|n><around*|[|k\<cdot\>G<around*|(|z<rsub|e>-z<rsub|+n>|)>+G<around*|(|z<rsub|e>+z<rsub|+n>|)>+G<around*|(|z<rsub|e>-z<rsub|-n>|)>+<frac|1|k>\<cdot\>G<around*|(|z<rsub|e>+z<rsub|-n>|)>|]>>>>>
  </eqnarray>

  <\eqnarray>
    <tformat|<table|<row|<cell|k>|<cell|=>|<cell|<frac|\<epsilon\><rsub|2>-\<epsilon\><rsub|3>|\<epsilon\><rsub|2>+\<epsilon\><rsub|3>>>>|<row|<cell|\<lambda\>>|<cell|=>|<cell|k\<cdot\><frac|\<epsilon\><rsub|2>-\<epsilon\><rsub|1>|\<epsilon\><rsub|2>+\<epsilon\><rsub|1>>>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|z<rsub|+n>>|<cell|=>|<cell|a<rsub|+n>>>|<row|<cell|>|<cell|=>|<cell|-z<rsub|h>+2*n
    w>>|<row|<cell|>|<cell|=>|<cell|-z<rsub|h>+4 n
    l>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|z<rsub|-n>>|<cell|=>|<cell|a<rsub|-n>>>|<row|<cell|>|<cell|=>|<cell|z<rsub|h>+2
    n w>>|<row|<cell|>|<cell|=>|<cell|z<rsub|h>+4 n l.>>>>
  </eqnarray>

  \;

  \;

  \;
</body>

<initial|<\collection>
</collection>>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|1.1|?>>
    <associate|auto-3|<tuple|1.2|?>>
    <associate|auto-4|<tuple|2|?>>
    <associate|auto-5|<tuple|3|?>>
    <associate|auto-6|<tuple|3.1|?>>
    <associate|auto-7|<tuple|3.2|?>>
    <associate|auto-8|<tuple|3.3|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Hamiltonian>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Kinetic terms
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>Potential terms
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Trial
      wave function> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Matrix
      elements> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Kinetic terms
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Potential terms
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|1tab>|3.3<space|2spc>One barrier case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>
    </associate>
  </collection>
</auxiliary>