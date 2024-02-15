#include <ibex.h>

using namespace std;
using namespace ibex;

int main(int argc, char* argv[]){

  // Dérivée symbolique directement sortie de la doc
  Function f("x","y","z","x*y*z");
	Function df(f,Function::DIFF);
	cout << "df=" << df << endl;

  // Il est aussi possible de dériver uniquement par rapport à certaines variables, par exemple ici, on va dériver par rapport à x et z
  Dim d(1,1); // Build the dimension of a scalar
  const ExprSymbol& x = ExprSymbol::new_("x",d); // nouvel argument de dimension (1,1)
  const ExprSymbol& y = ExprSymbol::new_("y",d);
  const ExprSymbol& z = ExprSymbol::new_("z",d);
  const Array<const ExprSymbol> new_args(x,y,z); // Définition des nouveaux arguments pour l'expression de la dérivée

  const Array<const ExprSymbol>& f_args(f.args()); // Arguments de f

  Array<const ExprSymbol> diff_args(2); // Arguments à dériver (x et z)
  diff_args.set_ref(0,f_args[0]); // Ajout du premier argument, ici f.args()[0] correspond à x
  diff_args.set_ref(1,f_args[2]); // Ajout du second argument, ici f.args()[2] correspond à z

  ExprDiff df_xz(f_args,new_args); // Définition de la classe ExprDiff
  const ExprNode& expr_df_xz = df_xz.diff(f.expr(),diff_args); // Calcul l'expression de la dérivée de f par rapport à x et z
  Function df_(new_args,expr_df_xz); // Création de la fonction dérivée
  cout << "df_=" << df_ << endl;

}
