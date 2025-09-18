# Contributing to **GEOppi**
We are thankful for contributions to **GEOppi**! Please refer to the information below for instructions on how to make contributions.

## Ways of contribution
- **Code adjustments and improvements**: New features, bug fixes, performance imporvements
- **Documentation**: README updates, code comments, tutorials
- **Feature requests**: Ask for new functionalities, e.g. on the [**GEOppi** Github discussion board](https://github.com/cvTHM/GEOppi/discussions)

## Contact
- **Issues**: Open an issue on the [**GEOppi** Github issue board](https://github.com/cvTHM/GEOppi/issues)
- **E-mail**: Contact the maintainer via constantin.voelzel@me.thm.de


## Making changes

### 1. Create a Feature Branch

```bash
# Always create a new branch for your changes
git checkout -b feature/descriptive-name
# or
git checkout -b fix/issue-number
```

### 2. Development Workflow

- **Pure Python changes**: Edit files directly and test with template files ([networkExtension_template.py](https://github.com/cvTHM/GEOppi/blob/main/networkExtension_template.py) and [networkModelling_template.py](https://github.com/cvTHM/GEOppi/blob/main/networkModelling_templatepy))
- **Documentation**: Update docstrings, README, or documentation files if appliccable

### 3. Commit Your Changes

```bash
# Stage your changes
git add .

# Commit with a descriptive message
git commit -m "feat: add new feature

- Implement EasterEgg feature
- Updated documentation
"
```

## Code Style Guidelines

### Python Code

- Use type hints where appropriate
- Add docstrings to all public functions and classes (sphinx, sphinx-notypes)

## Pull Request Process

### Before Submitting

- [ ] Code follows project style guidelines
- [ ] Documentation is updated (if applicable)
- [ ] Optional, if useful: New tests added for new functionality

### Pull Request Checklist

**Create Pull Request**
```bash
# Push your branch
git push origin feature/your-feature-name
```
   
**PR Description Should Include:**
- Clear description of changes
- Issue number (if applicable): `Fixes #123`
- Breaking changes (if any)

### Reviews
- A short-term review process can currently NOT be guaranteed.

## Reporting Issues

### Bug Reports

Use the [bug report template](https://github.com/cvTHM/GEOppi/issues/new?template=bug_report.md) and include:
- Python version and operating system
- Complete error traceback
- Minimal code example to reproduce
- Expected vs. actual behavior

## Suggesting Enhancements

Use the [enhancement template](https://github.com/cvTHM/GEOppi/issues/new?template=enhancement.md) and include:
- Clear description of the proposed feature
- Use case and motivation
- Possible implementation approach
- Performance considerations (if applicable)

### Security Issues

For security-related issues, email constantin.voelzel@me.thm.de directly instead of opening a public issue.