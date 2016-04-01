<TeXmacs|1.99.2>

<style|generic>

<\body>
  <doc-data|<doc-title|Electromagnetic fields in heavy-ion
  collisions>|<doc-author|<author-data|<author-name|Chun
  SHEN>|<\author-affiliation>
    McGill University
  </author-affiliation>>>>

  <abstract-data|<abstract|In this short note, we derived the space-time
  evolution of electromagnetic fields in relativistic heavy-ion collisions.>>

  Here, we derive the space-time evolution of the electric and magentic
  fields in relativistic heavy-ion colliisons originated from the fast moving
  charges inside the spectators and participants of the two colliding nuclei.

  In the local rest frame of a moving charge, at position
  <math|<around*|(|x<rsub|i><rprime|'>,y<rsub|i><rprime|'>|)>> with velocity
  <math|<math-bf|v>> in the lab frame, only electric field exists,

  <\equation>
    <math-bf|E><rprime|'><rsub|><around*|(|<math-bf|r><rprime|'>|)>=<frac|q|4\<pi\>><frac|<math-bf|r><rprime|'>|r<rprime|'><rsup|3>>,
  </equation>

  where <math|<math-bf|r><rprime|'>=<around*|(|x<rprime|'>-x<rsub|i><rprime|'>,y<rprime|'>-y<rsub|i><rprime|'>,z<rprime|'>|)>>.
  And

  <\equation>
    <math-bf|B><rprime|'><around*|(|<math-bf|r><rprime|'>|)>=0.
  </equation>

  In heavy-ion collisoins, <math|<math-bf|v> = v<rsub|z><math-bf|e<rsub|z>>>.
  We can boost the <math|<math-bf|E><rprime|'>> and
  <math|<math-bf|B><rprime|'>> fields back to lab frame,

  <\equation>
    <math-bf|E>=\<gamma\><around*|(|<math-bf|E><rprime|'>-<math-bf|\<beta\>>\<times\><math-bf|B><rprime|'>|)>-<frac|\<gamma\><rsup|2>|\<gamma\>+1><math-bf|\<beta\>><around*|(|<math-bf|\<beta\>>\<cdot\><math-bf|E><rprime|'>|)>=\<gamma\><math-bf|E><rprime|'>-<frac|\<gamma\><rsup|2>|\<gamma\>+1><math-bf|\<beta\>><around*|(|<math-bf|\<beta\>>\<cdot\><math-bf|E><rprime|'>|)>,
  </equation>

  <\equation>
    <math-bf|B> = \<gamma\><around*|(|<math-bf|B><rprime|'>+<math-bf|\<beta\>>\<times\><math-bf|E><rprime|'>|)>-<frac|\<gamma\><rsup|2>|\<gamma\>+1><math-bf|\<beta\>><around*|(|<math-bf|\<beta\>>\<cdot\><math-bf|B><rprime|'>|)>=\<gamma\><math-bf|\<beta\>>\<times\><math-bf|E><rprime|'>.
  </equation>

  with <math|<math-bf|\<beta\>>=<math-bf|v>=v<rsub|z><math-bf|e<rsub|z>>>.\ 

  \;

  From arXiv:1602.02223, the authors compute the EM fields in the presence of
  finite electric conductivity <math|\<sigma\>> and finite chiral magnetic
  conductivity <math|\<sigma\><rsub|\<chi\>>>. Assuming
  <math|\<sigma\><rsub|\<chi\>>\<ll\> \<sigma\>>, the authors were able to
  get analytic solution for the magnetic field. The authors futher assume the
  charge particle is highly relativistic, <math|\<gamma\>\<gg\>1>, then they
  can get analytic expression for the electric fields. The following is their
  results,

  <\equation>
    B<rsub|\<phi\>><around*|(|t,<math-bf|x>|)>=<frac|Q|4\<pi\>><frac|v\<gamma\>x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|(|1+<frac|\<sigma\>v\<gamma\>|2><sqrt|\<Delta\>>|)>e<rsup|A>
  </equation>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|B<rsub|r><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|-\<sigma\><rsub|\<chi\>><frac|Q|8\<pi\>><frac|v\<gamma\>x<rsub|\<perp\>>|\<Delta\><rsup|3/2>>
    <around*|[|v\<gamma\><around*|(|t-<frac|z|v>|)>+A<sqrt|\<Delta\>>|]>e<rsup|A>,<eq-number>>>|<row|<cell|B<rsub|z><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|\<sigma\><rsub|\<chi\>><frac|Q|8\<pi\>><frac|v\<gamma\>|\<Delta\><rsup|3/2>><around*|[|v<rsup|2>\<gamma\><rsup|2><around*|(|t-<frac|z|v>|)><rsup|2><around*|(|1+<frac|\<sigma\>v\<gamma\>|2><sqrt|\<Delta\>>|)>+\<Delta\><around*|(|1-<frac|\<sigma\>v\<gamma\>|2><sqrt|\<Delta\>>|)>|]>e<rsup|A>,<eq-number>>>>>
  </eqnarray*>

  where

  <\equation>
    \<Delta\>=v<rsup|2>\<gamma\><rsup|2><around*|(|t-z/v|)><rsup|2>+x<rsub|\<perp\>><rsup|2>
  </equation>

  and

  <\equation>
    A =<frac|\<sigma\>v\<gamma\>|2><around*|[|v\<gamma\><around*|(|t-<frac|z|v>|)>-<sqrt|\<Delta\>>|]>.
  </equation>

  And the expressions for the electric fields are,

  <\eqnarray*>
    <tformat|<table|<row|<cell|E<rsub|\<phi\>><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|\<sigma\><rsub|\<chi\>><frac|Q|8\<pi\>><frac|v<rsup|2>\<gamma\><rsup|2>x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|[|v\<gamma\><around*|(|t-<frac|z|v>|)>+A<sqrt|\<Delta\>>|]>e<rsup|A><eq-number>>>|<row|<cell|E<rsub|r><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|<frac|Q|4\<pi\>>e<rsup|A><around*|{|<frac|\<gamma\>x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|(|1+<frac|\<sigma\>v\<gamma\>|2><sqrt|\<Delta\>>|)>-<frac|\<sigma\>|v
    x<rsub|\<perp\>>><around*|[|1+<frac|v\<gamma\><around*|(|t-z/v|)>|<sqrt|\<Delta\>>>|]>|}><eq-number>>>|<row|<cell|E<rsub|z><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|<frac|Q|4\<pi\>><around*|{|e<rsup|-A><frac|1|\<Delta\><rsup|3/2>><around*|[|v\<gamma\><around*|(|t-<frac|z|v>|)>+A<sqrt|\<Delta\>>+<frac|\<sigma\>\<gamma\>|v>\<Delta\>|]>+<frac|\<sigma\><rsup|2>|v<rsup|2>>e<rsup|-\<sigma\><around*|(|t-z/v|)>>\<Gamma\><around*|(|0,-A|)>|}><eq-number>>>>>
  </eqnarray*>

  where <math|\<Gamma\><around*|(|0\<nocomma\>,-A|)>> is the incomplete gamma
  function defined as <math|\<Gamma\><around*|(|a,z|)>=<big|int><rsub|z><rsup|\<infty\>>d
  t t<rsup|a-1>exp<around*|(|-t|)>>.

  Now, we can write <math|v> in terms of rapidity,

  <\equation*>
    v=tanh<around*|(|Y<rsub|b>|)>,\<gamma\>=cosh<around*|(|Y<rsub|b>|)>
  </equation*>

  So

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<Delta\>>|<cell|=>|<cell|v<rsup|2>\<gamma\><rsup|2><around*|(|t-z/v|)><rsup|2>+x<rsub|\<perp\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|sinh<rsup|2><around*|(|Y<rsub|b>|)><around*|(|\<tau\>cosh<around*|(|\<eta\><rsub|s>|)>-<frac|\<tau\>sinh<around*|(|\<eta\><rsub|s>|)>|tanh<around*|(|Y<rsub|b>|)>>|)><rsup|2>+x<rsub|\<perp\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|\<tau\>sinh<around*|(|Y<rsub|b>-\<eta\><rsub|s>|)>|)><rsup|2>+x<rsub|\<perp\>><rsup|2>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|A>|<cell|=>|<cell|<frac|\<sigma\>v\<gamma\>|2><around*|[|v\<gamma\><around*|(|t-<frac|z|v>|)>-<sqrt|\<Delta\>>|]>>>|<row|<cell|>|<cell|=>|<cell|<frac|\<sigma\>|2>sinh<around*|(|Y<rsub|b>|)><around*|[|\<tau\>sinh<around*|(|Y<rsub|b>-\<eta\><rsub|s>|)>-<sqrt|\<Delta\>>|]>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|e B<rsub|\<phi\>><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|<frac|e<rsup|2>|4\<pi\>>Q<frac|v\<gamma\>x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|(|1+<frac|\<sigma\>v\<gamma\>|2><sqrt|\<Delta\>>|)>e<rsup|A>>>|<row|<cell|>|<cell|=>|<cell|\<alpha\><rsub|EM>Q
    sinh<around*|(|Y<rsub|b>|)><frac|x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|(|1+<frac|\<sigma\>|2>sinh<around*|(|Y<rsub|b>|)><sqrt|\<Delta\>>|)>e<rsup|A>>>|<row|<cell|e
    B<rsub|r><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|-\<sigma\><rsub|\<chi\>><frac|e<rsup|2>|8\<pi\>>Q<frac|v\<gamma\>
    x<rsub|\<perp\>>|\<Delta\><rsup|3/2>>
    <around*|[|v\<gamma\><around*|(|t-<frac|z|v>|)>+A<sqrt|\<Delta\>>|]>e<rsup|A>>>|<row|<cell|>|<cell|=>|<cell|-\<alpha\><rsub|EM>Q
    <frac|\<sigma\><rsub|\<chi\>>|2> sinh<around*|(|Y<rsub|b>|)><frac|x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|[|\<tau\>sinh<around*|(|Y<rsub|b>-\<eta\><rsub|s>|)>+A<sqrt|\<Delta\>>|]>e<rsup|A>>>|<row|<cell|e
    B<rsub|z><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|\<sigma\><rsub|\<chi\>><frac|e
    Q|8\<pi\>><frac|v\<gamma\>|\<Delta\><rsup|3/2>><around*|[|v<rsup|2>\<gamma\><rsup|2><around*|(|t-<frac|z|v>|)><rsup|2><around*|(|1+<frac|\<sigma\>v\<gamma\>|2><sqrt|\<Delta\>>|)>+\<Delta\><around*|(|1-<frac|\<sigma\>v\<gamma\>|2><sqrt|\<Delta\>>|)>|]>e<rsup|A>>>|<row|<cell|>|<cell|=>|<cell|\<alpha\><rsub|EM>Q
    <frac|\<sigma\><rsub|\<chi\>>|2> sinh<around*|(|Y<rsub|b>|)><frac|e<rsup|A>|\<Delta\><rsup|3/2>>>>|<row|<cell|>|<cell|>|<cell|<space|4em>\<times\><around*|[|<around*|(|\<tau\>sinh<around*|(|Y<rsub|b>-\<eta\><rsub|s>|)>|\<nobracket\>><rsup|2><around*|(|1+<frac|\<sigma\>|2>sinh<around*|(|Y<rsub|b>|)><sqrt|\<Delta\>>|)>+\<Delta\><around*|(|1-<frac|\<sigma\>|2>sinh<around*|(|Y<rsub|b>|)><sqrt|\<Delta\>>|)>|]>>>>>
  </eqnarray*>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|e E<rsub|\<phi\>><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|\<sigma\><rsub|\<chi\>><frac|e
    Q|8\<pi\>><frac|v<rsup|2>\<gamma\><rsup|2>x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|[|v\<gamma\><around*|(|t-<frac|z|v>|)>+A<sqrt|\<Delta\>>|]>e<rsup|A>>>|<row|<cell|>|<cell|=>|<cell|\<alpha\><rsub|EM>Q
    <frac|\<sigma\><rsub|\<chi\>>|2>sinh<rsup|2><around*|(|Y<rsub|b>|)><frac|x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|[|\<tau\>sinh<around*|(|Y<rsub|b>-\<eta\><rsub|s>|)>+A<sqrt|\<Delta\>>|]>e<rsup|A>>>|<row|<cell|e
    E<rsub|r><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|<frac|e
    Q|4\<pi\>>e<rsup|A><around*|{|<frac|\<gamma\>x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|(|1+<frac|\<sigma\>v\<gamma\>|2><sqrt|\<Delta\>>|)>-<frac|\<sigma\>|v
    x<rsub|\<perp\>>><around*|[|1+<frac|v\<gamma\><around*|(|t-z/v|)>|<sqrt|\<Delta\>>>|]>|}>>>|<row|<cell|>|<cell|=>|<cell|\<alpha\><rsub|EM>Q
    e<rsup|A><around*|{|cosh<around*|(|Y<rsub|b>|)><frac|x<rsub|\<perp\>>|\<Delta\><rsup|3/2>><around*|(|1+<frac|\<sigma\>|2>sinh<around*|(|Y<rsub|b>|)><sqrt|\<Delta\>>|)>-<frac|\<sigma\>|tanh<around*|(|Y<rsub|b>|)>x<rsub|\<perp\>>><around*|[|1+<frac|\<tau\>sinh<around*|(|Y<rsub|b>-\<eta\><rsub|s>|)>|<sqrt|\<Delta\>>>|]>|}>>>|<row|<cell|e
    E<rsub|z><around*|(|t,<math-bf|x>|)>>|<cell|=>|<cell|<frac|e
    Q|4\<pi\>><around*|{|e<rsup|-A><frac|1|\<Delta\><rsup|3/2>><around*|[|v\<gamma\><around*|(|t-<frac|z|v>|)>+A<sqrt|\<Delta\>>+<frac|\<sigma\>\<gamma\>|v>\<Delta\>|]>+<frac|\<sigma\><rsup|2>|v<rsup|2>>e<rsup|-\<sigma\><around*|(|t-z/v|)>>\<Gamma\><around*|(|0,-A|)>|}>>>|<row|<cell|>|<cell|=>|<cell|\<alpha\><rsub|EM>Q<around*|{|<frac|e<rsup|-A>|\<Delta\><rsup|3/2>><around*|[|\<tau\>sinh<around*|(|Y<rsub|b>-\<eta\><rsub|s>|)>+A<sqrt|\<Delta\>>+\<sigma\><frac|cosh<around*|(|Y<rsub|b>|)>|tanh<around*|(|Y<rsub|b>|)>>\<Delta\>|]>|\<nobracket\>>>>|<row|<cell|>|<cell|>|<cell|<space|4em>+<frac|\<sigma\><rsup|2>|tanh<rsup|2><around*|(|Y<rsub|b>|)>>exp<around*|[|-\<sigma\><around*|\<nobracket\>|<around*|(|<frac|\<tau\>sinh<around*|(|Y<rsub|b>-\<eta\><rsub|s>|)>|sinh<around*|(|Y<rsub|b>|)>>|)>|]>\<Gamma\><around*|(|0,-A|)>|}>>>>>
  </eqnarray*>
</body>