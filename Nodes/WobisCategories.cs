using Synera.Core.Implementation.UI;
using Synera.Core.UI;

public static class WobisCategories
{
    public static readonly ICategory Aida = new Category("Aida", 200, "Aida", "Aida", "Stress field creation for a solved fea model.");
}

public static class WobisSubcategories
{
    // If you extend your own category, create a nested static class with your category's name:
    public static readonly ICategory Field = new Category("Field", 201, "Field");
}
