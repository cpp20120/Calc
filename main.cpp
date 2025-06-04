import calc;
import std;

void display_error(const calc::ArithmeticError error) {
  switch (error) {
    case calc::ArithmeticError::Overflow:
      std::println("Error: Arithmetic overflow");
      break;
    case calc::ArithmeticError::Underflow:
      std::println("Error: Arithmetic underflow");
      break;
    case calc::ArithmeticError::InvalidInput:
      std::println("Error: Invalid number format");
      break;
    case calc::ArithmeticError::OutOfRange:
      std::println("Error: Number out of range");
      break;
    case calc::ArithmeticError::DivisionByZero:
      std::println("Error: Division by zero");
      break;
    case calc::ArithmeticError::InvalidExpression:
      std::println("Error: Invalid expression format");
      break;
    case calc::ArithmeticError::UnbalancedParentheses:
      std::println("Error: Unbalanced parentheses");
      break;
    case calc::ArithmeticError::UnknownOperator:
      std::println("Error: Unknown operator");
      break;
    default:
      std::println("Error: Unknown error");
  }
}

void display_result(
    const std::expected<double, calc::ArithmeticError>& result) {
  if (!result) {
    display_error(result.error());
    return;
  }

  double value = *result;
  if (value == std::floor(value) && value <= std::numeric_limits<int>::max() &&
      value >= std::numeric_limits<int>::min()) {
    std::println("= {}", static_cast<int>(value));
  } else {
    std::println("= {:.6}", value);
  }
}

int main() {
  std::println("Safe Calculator (type 'exit' to quit)");
  std::println("Supported operations: + - * / () !");
  std::println(
      "Supported math functions: sin cos tan cot sec csc asin acos atan acot "
      "asec acsc log ln sqrt ");
  std::println("Example: (2 + 3) * 4 - 10 / 2");
  std::println("         -5 + 3");
  std::println("         2 * -3");

  while (true) {
    std::print("> ");
    std::string input;
    std::getline(std::cin, input);

    input.erase(input.begin(),
                std::find_if(input.begin(), input.end(),
                             [](int ch) { return !std::isspace(ch); }));
    input.erase(std::find_if(input.rbegin(), input.rend(),
                             [](int ch) { return !std::isspace(ch); })
                    .base(),
                input.end());

    if (input.empty()) continue;
    if (input == "exit" || input == "quit") break;

    auto result = calc::calculate<double>(input);
    display_result(result);
  }

  std::println("Goodbye!");
  return 0;
}